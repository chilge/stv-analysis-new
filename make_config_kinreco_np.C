#include "ConfigMakerUtils.hh"
#include "ResponseMatrixMaker.hh"
#include "SliceBinning.hh"

constexpr int DUMMY_BLOCK_INDEX = -1;

void make_config_kinreco_np() {

  // Keys are the reco STV ntuple branch names of interest. Values
  // are vectors of bin edges
  std::map< std::string, std::vector<double> > kinreco_bin_edge_map = {
    { "kin_reco_pmu_diff", {-1.5,-0.85,-0.13,0.2,3} },
    { "kin_reco_costhetamu_diff", {-2,-0.28,0.18,0.8,2} }
  };

  // Used for converting from the variable names given in the bin edge
  // map to the ones used by the SliceBinning object
  std::map< std::string, std::string > var_name_map = {
    { "kin_reco_pmu_diff", "reco #Deltap_{#mu}" },
    { "kin_reco_costhetamu_diff", "reco #Deltacos#theta_{#mu}" },
  };

  // Keys are the same reco variable branch expressions as above. Values
  // are bool pairs indicating whether an underflow and overflow bin
  // should be produced.
  std::map< std::string, std::pair<bool,bool> > kinreco_under_overflow_map = {
    { "kin_reco_pmu_diff", { true, true } },
    { "kin_reco_costhetamu_diff", { false, false } }
  };

  // Set up an initially empty container to hold the slice definitions. We'll
  // populate it in parallel with defining the bins themselves.
  SliceBinning sb;

  // Set the variables to use when defining phase-space slices
  sb.slice_vars_ = {
    { "reco #Deltap_{#mu}", "GeV/c", "reco $\\Deltap_{\\mu}$", "GeV$/c$" },
    { "reco #Deltacos#theta_{#mu}", "", "reco $\\Delta\\cos\\theta_{\\mu}$", "" },
    { "reco bin number", "", "reco bin number", "" }
  };

  std::string selection =
  "sel_CCNp0pi "
  " && kin_reco_pmu_diff!=9999 && kin_reco_costhetamu_diff!=9999"
  " && kin_reco_enu>0 && kin_reco_enu<6"
  " && kin_reco_pmu>0 && kin_reco_pmu<6 "
  " && np>1";

  std::string signal_def =
  "mc_is_signal"
  " && mc_kin_reco_pmu_diff != 9999 && mc_kin_reco_costhetamu_diff!=9999"
  " && mc_kin_reco_enu > 0 && mc_kin_reco_enu < 6"
  " && mc_kin_reco_pmu > 0 && mc_kin_reco_pmu < 6"
  " && mc_np > 1";

  // By construction, MC event categories 5-11 contain all beam-correlated
  // backgrounds. This list is therefore comprehensive apart from cosmic
  // overlay stuff which is directly measured in a dedicated sample.
  std::vector< std::string > background_defs = {
    "category == 5", "category == 6", "category == 7", "category == 8",
    "category == 9", "category == 10", "category == 11"
  };

  std::vector< TrueBin > true_bins;
  std::vector< RecoBin > reco_bins;

  // Create separate blocks of bins for each kinematic variable
  int block_idx = -1;
  for ( const auto& pair : kinreco_bin_edge_map ) {
    // Start a new block of related bins
    ++block_idx;

    std::string reco_branchexpr = pair.first;
    std::string true_branchexpr = "mc_" + reco_branchexpr;

    std::string selection_tmp = selection;
    std::string signal_tmp = signal_def;
    if(pair.first == "kin_reco_pmu_diff"){
      selection_tmp += " && reco_muon_contained";
      signal_tmp += " && reco_muon_contained";
    }

    // Get the index for the "active" variable in the current block. We will
    // use it below to make a new slice while also defining the bins in the
    // block.
    const std::string& act_var_name = var_name_map.at( reco_branchexpr );
    int act_var_idx = find_slice_var_index( act_var_name, sb.slice_vars_ );

    const auto flag_pair = kinreco_under_overflow_map.at( reco_branchexpr );
    bool needs_underflow_bin = flag_pair.first;
    bool needs_overflow_bin = flag_pair.second;

    // Require at least two bin edges to be present in the input vector.
    // Any variables for which this is not true will be skipped entirely.
    const auto& bin_edges = pair.second;

    size_t num_edges = bin_edges.size();
    size_t num_bins = 0u;
    if ( num_edges >= 2u ) num_bins = num_edges - 1u;
    else continue;

    // Before defining each bin, make a new Slice object and set up the
    // corresponding ROOT histogram within it
    auto& cur_slice = add_slice( sb, bin_edges, act_var_idx );

    // If needed, then create the underflow bin in both true and reco space
    if ( needs_underflow_bin ) {
      double var_underflow_max = bin_edges.front();

      std::stringstream true_ss;
      true_ss << signal_tmp << " && " << true_branchexpr
        << " < " << var_underflow_max;

      std::string true_bin_def = true_ss.str();
      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_idx );

      std::stringstream reco_ss;
      reco_ss << selection_tmp << " && " << reco_branchexpr
        << " < " << var_underflow_max;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry. Note that the ROOT histogram bin
      // indices are one-based, so the underflow bin is always at index zero.
      cur_slice.bin_map_[ 0 ].insert( ana_bin_idx );

      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_idx );
    }

    // Create the ordinary signal bins using the requested edges
    for ( size_t b = 0u; b < num_bins; ++b ) {

      double var_low = bin_edges.at( b );
      double var_high = bin_edges.at( b + 1u );

      std::stringstream true_ss;
      true_ss << signal_tmp
        << " && " << true_branchexpr << " >= " << var_low
        << " && " << true_branchexpr << " < "  << var_high;

      std::string true_bin_def = true_ss.str();

      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_idx );

      std::stringstream reco_ss;
      reco_ss << selection_tmp
        << " && " << reco_branchexpr << " >= " << var_low
        << " && " << reco_branchexpr << " < "  << var_high;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry. Note that the ROOT histogram bin
      // indices are one-based, so we correct for that in the line below.
      cur_slice.bin_map_[ b + 1 ].insert( ana_bin_idx );

      // Add the completed reco bin definition to the vector
      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_idx );

    } // loop over ordinary bins for the current variable

    // If needed, then create the overflow bin in both true and reco space
    if ( needs_overflow_bin ) {
      double var_overflow_min = bin_edges.back();

      std::stringstream true_ss;
      true_ss << signal_tmp << " && " << true_branchexpr
        << " >= " << var_overflow_min;

      std::string true_bin_def = true_ss.str();
      true_bins.emplace_back( true_bin_def, kSignalTrueBin, block_idx );

      std::stringstream reco_ss;
      reco_ss << selection_tmp << " && " << reco_branchexpr
        << " >= " << var_overflow_min;

      std::string reco_bin_def = reco_ss.str();

      // Here we use a trick: the current analysis bin index is equal
      // to the size of the reco_bins vector before we add the new element.
      size_t ana_bin_idx = reco_bins.size();
      // Here's another trick: the call to operator[]() below will create
      // a new map entry if needed. We then insert the current analysis
      // bin index into the map entry. Note that the ROOT histogram bin
      // indices are one-based, so the overflow bin has an index equal to
      // the number of bin edges.
      cur_slice.bin_map_[ bin_edges.size() ].insert( ana_bin_idx );

      // Add the completed reco bin definition to the vector
      reco_bins.emplace_back( reco_bin_def, kOrdinaryRecoBin, block_idx );
    }

  } // loop over kinematic variables

  // Add a single set of true bins for the background categories of interest
  for ( const auto& bdef : background_defs ) {
    true_bins.emplace_back( bdef, kBackgroundTrueBin, DUMMY_BLOCK_INDEX );
  }

  // Create a slice showing all blocks together as a function of bin number
  int num_reco_bins = reco_bins.size();
  int bin_number_var_idx = find_slice_var_index( "reco bin number",
    sb.slice_vars_ );

  auto& bin_num_slice = add_slice( sb, num_reco_bins, 0, num_reco_bins,
    bin_number_var_idx );
  for ( int ab = 0; ab < num_reco_bins; ++ab ) {
    // The ROOT histogram bins are one-based, so we correct for this here
    bin_num_slice.bin_map_[ ab + 1 ].insert( ab );
  }

  // Dump this information to the output file
  std::ofstream out_file( "myconfig_kinreco_np.txt" );
  out_file << "kinreco_all" << '\n';
  out_file << "stv_tree\n";
  out_file << true_bins.size() << '\n';
  for ( const auto& tb : true_bins ) out_file << tb << '\n';

  out_file << reco_bins.size() << '\n';
  for ( const auto& rb : reco_bins ) out_file << rb << '\n';

  // Also write a SliceBinning configuration file
  std::ofstream sb_file( "mybins_kinreco_np.txt" );
  sb_file << sb;
}
