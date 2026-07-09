// SystematicsCalculator caches its POT-summed universes the first time it
// processes a given univmake output file, in a 'total_<mangled config path>'
// TDirectoryFile written alongside the per-file universe histograms (see
// SystematicsCalculator::SystematicsCalculator in
// src/utils/SystematicsCalculator.cxx). On every later run, it only checks
// whether that directory *exists* -- never whether file_properties.txt, the
// systematics config, or the selection code have changed since -- so a cache
// built before you finished editing those inputs will silently keep serving
// stale (sometimes all-zero) numbers forever.
//
// This macro removes every 'total_*' key (at any depth) from the univmake
// output file, forcing the next tool that opens the file (SlicePlots,
// Unfolder, ...) to rebuild it from the per-file histograms.
//
// Note this does NOT use TDirectory::Delete()/rmdir(). Verified empirically
// against both toy files and a real univmake output file: once a
// TDirectoryFile key has had objects written into it (which is always true
// for these caches), neither call actually removes the key from a file
// opened in "update" mode -- they report success, but the key is still
// found by a later GetObject() call. The only reliable way to remove data
// from a ROOT file is to copy everything you want to keep into a new file
// and swap it in for the original, which is what this macro does.
//
// Usage: root -l -b -q 'invalidate_total_cache.C("univmake_output.root")'

namespace {

  // Cheap read-only scan: is there anything to do? Avoids paying for a full
  // file copy on the (common) case where the cache is already clean.
  bool has_stale_total_dir( const std::string& univ_file ) {
    TFile f( univ_file.c_str(), "read" );
    if ( f.IsZombie() ) {
      std::cerr << "invalidate_total_cache: could not open " << univ_file
        << '\n';
      return false;
    }

    TIter next_key( f.GetListOfKeys() );
    TKey* key = nullptr;
    while ( ( key = (TKey*)next_key() ) ) {
      if ( std::string(key->GetClassName()) != "TDirectoryFile" ) continue;

      TDirectoryFile* top_dir = nullptr;
      f.GetObject( key->GetName(), top_dir );
      if ( !top_dir ) continue;

      TIter next_inner( top_dir->GetListOfKeys() );
      TKey* inner_key = nullptr;
      while ( ( inner_key = (TKey*)next_inner() ) ) {
        std::string name = inner_key->GetName();
        if ( name.rfind( "total_", 0 ) == 0 ) return true;
      }
    }
    return false;
  }

  // Recursively copies 'source' into 'dest', dropping any key (at any
  // depth) whose name starts with 'total_'.
  void copy_dir_skip_stale_totals( TDirectory* source, TDirectory* dest ) {

    TIter next_key( source->GetListOfKeys() );
    TKey* key = nullptr;
    while ( ( key = (TKey*)next_key() ) ) {

      std::string name = key->GetName();
      if ( name.rfind( "total_", 0 ) == 0 ) {
        std::cout << "invalidate_total_cache: dropping " << name
          << " while rebuilding the file\n";
        continue;
      }

      TClass* cl = TClass::GetClass( key->GetClassName() );
      if ( !cl ) continue;

      if ( cl->InheritsFrom( TDirectory::Class() ) ) {
        TDirectory* source_sub = (TDirectory*)key->ReadObj();
        dest->cd();
        TDirectory* dest_sub = dest->mkdir( name.c_str() );
        copy_dir_skip_stale_totals( source_sub, dest_sub );
        delete source_sub;
      }
      else if ( cl->InheritsFrom( TObject::Class() ) ) {
        TObject* obj = key->ReadObj();
        dest->cd();
        obj->Write();
        delete obj;
      }
      else {
        // Non-TObject payloads (e.g. the raw std::string keys written by
        // UniverseMaker via TDirectory::WriteObject -- ntuple_name,
        // true_bin_spec, reco_bin_spec, sel_for_categ). TKey::ReadObj()
        // only handles TObject-derived classes and silently returns
        // nullptr for these, which previously caused a segfault here.
        // ReadObjectAny()/WriteObjectAny() are the generic equivalents.
        void* obj = key->ReadObjectAny( cl );
        dest->cd();
        dest->WriteObjectAny( obj, cl, name.c_str() );
        cl->Destructor( obj );
      }
    } // keys

    dest->SaveSelf( true );
  }

} // anonymous namespace

void invalidate_total_cache( std::string univ_file ) {

  if ( !has_stale_total_dir( univ_file ) ) {
    std::cout << "invalidate_total_cache: no cached 'total_*' directory"
      " found in " << univ_file << " -- nothing to do.\n";
    return;
  }

  std::string tmp_file = univ_file + ".invalidate_tmp";

  {
    TFile f_in( univ_file.c_str(), "read" );
    TFile f_out( tmp_file.c_str(), "recreate" );

    copy_dir_skip_stale_totals( &f_in, &f_out );

    f_out.Close();
    f_in.Close();
  }

  if ( gSystem->Rename( tmp_file.c_str(), univ_file.c_str() ) != 0 ) {
    std::cerr << "invalidate_total_cache: failed to replace " << univ_file
      << " with the rebuilt copy at " << tmp_file
      << " -- please check disk space / permissions and retry\n";
  }
}
