package com.compomics.cli.identification_parameters;

import com.compomics.software.cli.CommandLineUtils;
import com.compomics.util.experiment.identification.modification.ModificationLocalizationScore;
import com.compomics.util.experiment.identification.spectrum_annotation.AnnotationParameters;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import org.apache.commons.cli.Options;

/**
 * Enum class specifying the SearchParameter command line option parameters to
 * create a SearchParameters object
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public enum IdentificationParametersCLIParams {

    //////////////////////////////////////////////////////////////////////////////////////////
    // IMPORTANT: Any change here must be reported in the wiki: 
    // https://github.com/compomics/compomics-utilities/wiki/IdentificationParametersCLI
    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////
    // General parameters
    //////////////////////////////////
    DESCRIPTION("description", "The description of the search parameters, short summary by default.", false, true),
    USAGE("usage", "Lists the available options.", false, false),
    USAGE_2("h", "Lists the available options.", false, false),
    USAGE_3("help", "Lists the available options.", false, false),
    MODS("mods", "Lists the available modifications.", false, false),
    ENZYMES("enzymes", "Lists the available enzymes.", false, false),
    IDENTIFICATION_PARAMETERS("id_params", "The identification parameters file to use.", false, true),
    OUT("out", "The destination Identification Parameters file (.par).", true, true),
    //////////////////////////////////
    // Search Parameters
    //////////////////////////////////
    PREC_PPM("prec_ppm", "Precursor ion tolerance unit: ppm (1) or Da (0), default is '1'.", false, true),
    FRAG_PPM("frag_ppm", "Fragment ion tolerance unit: ppm (1) or Da (0), default is '0'.", false, true),
    PREC_TOL("prec_tol", "Precursor ion mass tolerance, default is '10'.", false, true),
    FRAG_TOL("frag_tol", "Fragment ion mass tolerance, default is '0.5'.", false, true),
    DIGESTION("digestion", "The type of digestion to consider: " + DigestionParameters.CleavageParameter.getCommandLineDescription() + ". Default is 0.", false, true),
    ENZYME("enzyme", "Enzyme used, default is 'Trypsin'. If more than one enzyme was used, please provide them as comma separated list with quotes, e.g. \"Trypsin, Glu-C\". See EnzymesCLI to list and edit the enzymes.", false, true),
    MC("mc", "Number of allowed missed cleavages, default is '2'. If more than one enzyme was used, please provide the missed cleavages for every enzyme in the same order as comma separated list with quotes, e.g. \"2, 1\".", false, true),
    SPECIFICITY("specificity", "Specificity of the enzyme." + DigestionParameters.Specificity.getCommandLineDescription() + ". If more than one enzyme was used, please provide the missed cleavages for every enzyme in the same order as comma separated list with quotes, e.g. \"0, 1\".", false, true),
    FIXED_MODS("fixed_mods", "Fixed modifications as comma separated list, e.g., \"Oxidation of M, Phosphorylation of S\"", false, true),
    VARIABLE_MODS("variable_mods", "Variable modifications as comma separated list, e.g., \"Oxidation of M, Phosphorylation of S\". See ModificationsCLI to list and edit the enzymes", false, true),
    MIN_CHARGE("min_charge", "Minimal charge to search for, default is '2'.", false, true),
    MAX_CHARGE("max_charge", "Maximal charge to search for, default is '4'.", false, true),
    FI("fi", "Type of forward ion searched, default is 'b'. If more than one ion should be used, please provide them as comma separated list with quotes, e.g. \"a, b\".", false, true),
    RI("ri", "Type of rewind ion searched, default is 'y'. If more than one ion should be used, please provide them as comma separated list with quotes, e.g. \"y, z\".", false, true),
    MIN_ISOTOPE("min_isotope", "Minimal precursor isotope, default is '0'.", false, true),
    MAX_ISOTOPE("max_isotope", "Maximal precursor isotope, default is '1'.", false, true),
    //////////////////////////////////
    // X!Tandem specific parameters
    //////////////////////////////////
    XTANDEM_DYNAMIC_RANGE("xtandem_dynamic_range", "X!Tandem 'spectrum, dynamic range' option, default is '100'.", false, true),
    XTANDEM_NPEAKS("xtandem_npeaks", "X!Tandem 'spectrum, total peaks' option, default is '50'.", false, true),
    XTANDEM_MIN_FRAG_MZ("xtandem_min_frag_mz", "X!Tandem 'spectrum, minimum fragment mz' option, default is '200'.", false, true),
    XTANDEM_MIN_PEAKS("xtandem_min_peaks", "X!Tandem 'spectrum, minimum peaks' option, default is '5'.", false, true),
    XTANDEM_NOISE_SUPPRESSION("xtandem_noise_suppr", "X!Tandem 'spectrum, use noise suppression' option. 1: true, 0: false, default is '0'.", false, true),
    XTANDEM_MIN_PREC_MASS("xtandem_min_prec_mass", "X!Tandem 'spectrum, minimum parent m+h' option, default is '500'.", false, true),
    XTANDEM_PARENT_MONOISOTOPIC_MASS_ISOTOPE_ERROR("xtandem_parent_isotope_correction", "X!Tandem 'spectrum, parent monoisotopic mass isotope error. 1: true, 0: false, default is '1'.", false, true),
    XTANDEM_QUICK_ACETYL("xtandem_quick_acetyl", "X!Tandem 'protein, quick acetyl' option. 1: true, 0: false, default is '1'.", false, true),
    XTANDEM_QUICK_PYRO("xtandem_quick_pyro", "X!Tandem 'protein, quick pyrolidone' option. 1: true, 0: false, default is '1'.", false, true),
    XTANDEM_STP_BIAS("xtandem_stp_bias", "X!Tandem 'protein, stP bias' option. 1: true, 0: false, default is '0'.", false, true),
    XTANDEM_PTM_COMPLEXITY("xtandem_ptm_complexity", "X!Tandem 'protein, ptm complexity' option, default is '6.0'.", false, true),
    XTANDEM_REFINE("xtandem_refine", "X!Tandem 'refine' option. 1: true, 0: false, default is '1'.", false, true),
    XTANDEM_REFINE_EVALUE("xtandem_refine_evalue", "X!Tandem 'refine, maximum valid expectation value' option, default is '0.01'.", false, true),
    XTANDEM_REFINE_UNANTICIPATED_CLEAVAGE("xtandem_refine_unc", "X!Tandem 'refine, unanticipated cleavage' option. 1: true, 0: false, default is '1'.", false, true),
    XTANDEM_REFINE_SEMI("xtandem_refine_semi", "X!Tandem 'refine, cleavage semi' option. 1: true, 0: false, default is '0'.", false, true),
    XTANDEM_REFINE_POTENTIAL_MOD_FULL_REFINEMENT("xtandem_refine_pot", "X!Tandem 'refine, use potential modifications for full refinement' option. 1: true, 0: false, default is '0'.", false, true),
    XTANDEM_REFINE_POINT_MUTATIONS("xtandem_refine_p_mut", "X!Tandem 'refine, point mutations' option. 1: true, 0: false, default is '0'.", false, true),
    XTANDEM_REFINE_SNAPS("xtandem_refine_snaps", "X!Tandem 'refine, sNAps' option. 1: true, 0: false, default is '1'.", false, true),
    XTANDEM_REFINE_SPECTRUM_SYNTHESIS("xtandem_refine_spec_synt", "X!Tandem 'refine, spectrum synthesis' option. 1: true, 0: false, default is '1'.", false, true),
    XTANDEM_EVALUE("xtandem_evalue", "X!Tandem 'output, maximum valid expectation value' option, default is '0.01'.", false, true),
    XTANDEM_OUTPUT_RESULTS("xtandem_output_results", "X!Tandem 'output, results' option (all|valid|stochastic), default is 'all'.", false, true),
    XTANDEM_OUTPUT_PROTEINS("xtandem_output_proteins", "X!Tandem 'output, proteins' option. 1: true, 0: false, default is '1'.", false, true),
    XTANDEM_OUTPUT_SEQUENCES("xtandem_output_sequences", "X!Tandem 'output, sequences' option. 1: true, 0: false, default is '0'.", false, true),
    XTANDEM_OUTPUT_SPECTRA("xtandem_output_spectra", "X!Tandem 'output, spectra' option. 1: true, 0: false, default is '1'.", false, true),
    XTANDEM_OUTPUT_HISTOGRAMS("xtandem_output_histograms", "X!Tandem 'output, histograms' option. 1: true, 0: false, default is '0'.", false, true),
    XTANDEM_SKYLINE("xtandem_skyline_path", "X!Tandem 'spectrum, skyline path' option.", false, true),
    //////////////////////////////////
    // MS-GF+ specific parameters
    //////////////////////////////////
    MSGF_DECOY("msgf_decoy", "MS-GF+ search decoys option, 1: true, 0: false, default is '0'.", false, true),
    MSGF_INSTRUMENT("msgf_instrument", "MS-GF+ instrument id option, 0: Low-res LCQ/LTQ, 1: Orbitrap/FTICR, 2: TOF, 3: Q-Exactive (Default).", false, true),
    MSGF_FRAGMENTATION("msgf_fragmentation", "MS-GF+ fragmentation id option, 0: Automatic: as written in the spectrum or CID if no info, 1: CID, 2: ETD, 3: HCD (Default). 4: UVPD.", false, true),
    MSGF_PROTOCOL("msgf_protocol", "MS-GF+ protocol id option, 0: Automatic (Default, true), 1: Phosphorylation, 2: iTRAQ, 3: iTRAQPhospho, 4: TMT, 5: Standard.", false, true),
    MSGF_TERMINI("msgf_termini", "MS-GF+ number of tolerable termini, e.g. 0: Non required, 1: At least one, 2: Both, default is '2'.", false, true),
    MSGF_MIN_PEP_LENGTH("msgf_min_pep_length", "MS-GF+ minumum peptide length, default is '8'.", false, true),
    MSGF_MAX_PEP_LENGTH("msgf_max_pep_length", "MS-GF+ maximum peptide length, default is '30'.", false, true),
    MSGF_PTMS("msgf_num_ptms", "MS-GF+ max number of PTMs per peptide, default is '2'.", false, true),
    MSGF_NUM_MATCHES("msgf_num_matches", "MS-GF+ maximum number of spectrum matches, default is '10'.", false, true), // @TODO: find an optimal default
    MSGF_ADDITIONAL("msgf_additional", "MS-GF+ additional output, 0: output basic scores only (Default, true), 1: output additional features.", false, true),
    MSGF_TASKS("msgf_num_tasks", "MS-GF+ number of tasks as an integer, default: internally calculated based on inputs", false, true),
    //////////////////////////////////
    // MS Amanda specific parameters
    //////////////////////////////////
    MS_AMANDA_DECOY("ms_amanda_decoy", "MS Amanda generate decoys option, 0: false, 1: true, default is '0'.", false, true),
    MS_AMANDA_DECOY_RANKING("ms_amanda_decoy_ranking", "MS Amanda decoys ranking options, 0: shared rank, 1: rank target and decoy separately, default is '1'.", false, true),
    MS_AMANDA_INSTRUMENT("ms_amanda_instrument", "MS Amanda instrument id option. Available enzymes are listed in the GUI. (Note: case sensitive.).", false, true),
    MS_AMANDA_MAX_RANK("ms_amanda_max_rank", "MS Amanda maximum rank, default is '10'.", false, true),
    MS_AMANDA_MONOISOTOPIC("ms_amanda_mono", "MS Amanda use monoisotopic mass values, 0: false, 1: true, default is '1'.", false, true),
    MS_AMANDA_PERFORM_DEISOTOPING("ms_amanda_perform_deisotoping", "MS Amanda perform deisotoping, 0: false, 1: true, default is '1'.", false, true),
    MS_AMANDA_MAX_MOD("ms_amanda_max_mod", "MS Amanda maximum number of occurrences of a specific modification on a peptide (0-10), default is '3'", false, true),
    MS_AMANDA_MAX_VAR_MOD("ms_amanda_max_var_mod", "MS Amanda maximum number of variable modifications per peptide (0-10), default is '4'", false, true),
    MS_AMANDA_MAX_MOD_SITES("ms_amanda_max_mod_sites", "MS Amanda maximum number of potential modification sites per modification per peptide (0-20), default is '6'", false, true),
    MS_AMANDA_MAX_NEUTRAL_LOSSES("ms_amanda_max_neutraul_losses", "MS Amanda maximum number of water and ammonia losses per peptide (0-5), default is '1'", false, true),
    MS_AMANDA_MAX_NEUTRAL_LOSSES_MODIFICATIONS("ms_amanda_max_neutral_losses_mod", "MS Amanda maximum number identical modification specific losses per peptide (0-5), default is '2'", false, true),
    MS_AMANDA_MIN_PEPTIDE_LENGTH("ms_amanda_min_pep_length", "MS Amanda minimum peptide length, default is '6'", false, true),
    MS_AMANDA_MAX_PEPTIDE_LENGTH("ms_amanda_max_pep_length", "MS Amanda maximum peptide length, default is '30'", false, true),
    MS_AMANDA_LOADED_PROTEINS("ms_amanda_loaded_proteins", "MS Amanda maximum number of proteins loaded into memory (1000-500000), default is '100000'", false, true),
    MS_AMANDA_LOADED_SPECTRA("ms_amanda_loaded_spectra", "MS Amanda maximum number of spectra loaded into memory (1000-500000), default is '2000'", false, true),
    MS_AMANDA_OUTPUT_FORMAT("ms_amanda_output", "MS Amanda output format option, csv or mzIdentML, default is 'csv'.", false, true),
    //////////////////////////////////
    // MyriMatch specific parameters
    //////////////////////////////////
    MYRIMATCH_MIN_PEP_LENGTH("myrimatch_min_pep_length", "MyriMatch minumum peptide length, default is '6'.", false, true),
    MYRIMATCH_MAX_PEP_LENGTH("myrimatch_max_pep_length", "MyriMatch maximum peptide length, default is '30'.", false, true),
    MYRIMATCH_MIN_PREC_MASS("myrimatch_min_prec_mass", "MyriMatch minumum precursor mass, default is '0.0'.", false, true),
    MYRIMATCH_MAX_PREC_MASS("myrimatch_max_prec_mass", "MyriMatch maximum precursor mass, default is '10000.0'.", false, true),
    MYRIMATCH_NUM_MATCHES("myrimatch_num_matches", "MyriMatch maximum number of spectrum matches, default is '10'.", false, true),
    MYRIMATCH_PTMS("myrimatch_num_ptms", "MyriMatch max number of PTMs per peptide, default is '2'.", false, true),
    MYRIMATCH_FRAGMENTATION("myrimatch_fragmentation", "MyriMatch fragmentation method, cid (b, y, true), etd (c, z*) or manual (a comma-separated list of [abcxyz] or z* (z+1, true), e.g. manual:b,y,z).", false, true),
    MYRIMATCH_TERMINI("myrimatch_termini", "MyriMatch number of enzymatic termini, e.g. 0: non-tryptic, 1: semi-tryptic, 2: fully-tryptic, default is '2'.", false, true),
    MYRIMATCH_SMART_PLUS_THREE("myrimatch_plus_three", "MyriMatch smart plus three option, 1: true, 0: false, default is '1'.", false, true),
    MYRIMATCH_XCORR("myrimatch_xcorr", "MyriMatch xcorr option, 1: true, 0: false, default is '0'.", false, true),
    MYRIMATCH_TIC_CUTOFF("myrimatch_tic_cutoff", "MyriMatch TIC cutoff percentage, default is '0.98'.", false, true),
    MYRIMATCH_INTENSTITY_CLASSES("myrimatch_intensity_classes", "MyriMatch number of intensity classes, default is '3'.", false, true),
    MYRIMATCH_CLASS_MULTIPLIER("myrimatch_class_multiplier", "MyriMatch class multiplier option, default is '2'.", false, true),
    MYRIMATCH_NUM_BATCHES("myrimatch_num_batches", "MyriMatch number of batches option, default is '50'.", false, true),
    MYRIMATCH_MAX_PEAK_COUNT("myrimatch_max_peak", "MyriMatch max number of peaks option, default is '100'.", false, true),
    MYRIMATCH_OUTPUT_FORMAT("myrimatch_output", "MyriMatch output format option, mzIdentML or pepXML, default is 'mzIdentML'.", false, true),
    //////////////////////////////////
    // OMSSA specific parameters
    //////////////////////////////////
    OMSSA_LOW_INTENSITY("omssa_low_intensity", "OMSSA low intensity cutoff as percentage of the most intense peak, default is '0.0'.", false, true),
    OMSSA_HIGH_INTENSITY("omssa_high_intensity", "OMSSA high intensity cutoff as percentage of the most intense peak, default is '0.2'.", false, true),
    OMSSA_INTENSITY_INCREMENT("omssa_intensity_incr", "OMSSA intensity cutoff increment, default is '0.0005'.", false, true),
    OMSSA_MIN_PEAKS("omssa_min_peaks", "OMSSA minimum number of peaks, integer, default is '4'.", false, true),
    OMSSA_REMOVE_PREC("omssa_remove_prec", "OMSSA eliminate charge reduced precursors in spectra option, 1: true, 0: false, default is '0'.", false, true),
    OMSSA_ESTIMATE_CHARGE("omssa_estimate_charge", "OMSSA estimate precursor charge option, 1: Use range, 0: Believe input file, default is '1'.", false, true),
    OMSSA_PLUS_ONE("omssa_plus_one", "OMSSA estimate plus one charge algorithmically option, 1: true, 0: false, default is '1'.", false, true),
    OMSSA_MAX_FRACTION("omssa_fraction", "OMSSA fraction of precursor m/z for charge charge one estimation, default is '0.95'.", false, true),
    OMSSA_PREC_PER_SPECTRUM("omssa_prec_per_spectrum", "OMSSA minimum number of precursors per spectrum, integer, default is '1'.", false, true),
    OMSSA_SCALE_PREC("omssa_scale_prec", "OMSSA scale precursor mass option, 1: true, 0: false, default is '1'.", false, true),
    OMSSA_SEQUENCES_IN_MEMORY("omssa_memory", "OMSSA map sequences in memory option, 1: true, 0: false, default is '1'.", false, true),
    OMSSA_METHIONINE("omssa_methionine", "OMSSA N-terminal methionine cleavage option, 1: true, 0: false, default is '1'.", false, true),
    // TODO: MINIMUM PRECURSOR CHARGE FOR MULTIPLY CHARGED FRAGMENTS
    OMSSA_NEUTRON("omssa_neutron", "Mass after which OMSSA should consider neutron exact mass, default is '1446.94'.", false, true),
    OMSSA_SINGLE_WINDOW_WIDTH("omssa_single_window_wd", "OMSSA single charge window width in Da, integer, default is '27'.", false, true),
    OMSSA_DOUBLE_WINDOW_WIDTH("omssa_double_window_wd", "OMSSA double charge window width in Da, integer, default is '14'.", false, true),
    OMSSA_SINGLE_WINDOW_PEAKS("omssa_single_window_pk", "OMSSA single charge window number of peaks, integer, default is '2'.", false, true),
    OMSSA_DOUBLE_WINDOW_PEAKS("omssa_double_window_pk", "OMSSA double charge window number of peaks, integer, default is '2'.", false, true),
    OMSSA_MIN_ANNOTATED_INTENSE_PEAKS("omssa_min_ann_int_pks", "OMSSA minimum number of annotated peaks among the most intense ones, integer, default is '6'.", false, true),
    OMSSA_MIN_ANNOTATED_PEAKS("omssa_min_annotated_peaks", "OMSSA minimum number of annotated peaks, integer, default is '2'.", false, true),
    OMSSA_MAX_LADDERS("omssa_max_ladders", "OMSSA maximum number of m/z ladders, integer, default is '128'.", false, true),
    OMSSA_MAX_FRAG_CHARGE("omssa_max_frag_charge", "OMSSA maximum fragment charge, integer, default is '2'.", false, true),
    OMSSA_POSITIVE_IONS("omssa_charge", "OMSSA fragment charge option, 1: plus, 0: minus, default is '1'.", false, true),
    OMSSA_FORWARD_IONS("omssa_forward", "OMSSA include first forward ion (b1) in search, 1: true, 0: false, default is '0'.", false, true),
    OMSSA_REWIND_IONS("omssa_rewind", "OMSSA search rewind (C-terminal) ions option, 1: true, 0: false, default is '1'.", false, true),
    OMSSA_MAX_FRAG_SERIES("omssa_max_frag_series", "OMSSA maximum fragment per series option, integer, default is '100'.", false, true),
    OMSSA_CORRELATION_CORRECTION("omssa_corr", "OMSSA use correlation correction score option, 1: true, 0: false, default is '1'.", false, true),
    OMSSA_CONSECUTIVE_ION_PROBABILITY("omssa_consecutive_p", "OMSSA consecutive ion probability, default is '0.5'.", false, true),
    OMSSA_HITLIST_LENGTH_CHARGE("omssa_hitlist_charge", "OMSSA number of hits per spectrum per charge, default is '30'.", false, true),
    OMSSA_ITERATIVE_SEQUENCE_EVALUE("omssa_it_sequence_evalue", "OMSSA e-value cutoff to consider a sequence in the iterative search 0.0 means all, default is '0.0'.", false, true),
    OMSSA_ITERATIVE_SPECTRUM_EVALUE("omssa_it_spectrum_evalue", "OMSSA e-value cutoff to consider a spectrum in the iterative search 0.0 means all, default is '0.01'.", false, true),
    OMSSA_ITERATIVE_REPLACE_EVALUE("omssa_it_replace_evalue", "OMSSA e-value cutoff to replace a hit in the iterative search 0.0 means keep best, default is '0.0'.", false, true),
    OMSSA_MIN_PEP_LENGTH("omssa_min_pep_length", "OMSSA minumum peptide length (semi-tryptic or no enzyme searches only).", false, true),
    OMSSA_MAX_PEP_LENGTH("omssa_max_pep_length", "OMSSA maximum peptide length (OMSSA semi-tryptic or no enzyme searches only).", false, true),
    OMSSA_MAX_EVALUE("omssa_max_evalue", "OMSSA maximal evalue considered, default is '100'.", false, true),
    OMSSA_HITLIST_LENGTH("omssa_hitlist_length", "OMSSA hitlist length, 0 means all, default is '10'.", false, true),
    OMSSA_FORMAT("omssa_format", "OMSSA output format. 0: omx, 1: csv, 2: pepXML, default is 'omx'.", false, true),
    //////////////////////////////////
    // Comet specific parameters
    //////////////////////////////////
    COMET_MIN_PEAKS("comet_min_peaks", "Comet min number of peaks for a spectrum, default is '10'.", false, true),
    COMET_MIN_PEAK_INTENSITY("comet_min_peak_int", "Comet minimal peak intensity, default is '0.0'.", false, true),
    COMET_REMOVE_PRECURSOR("comet_remove_prec", "Comet remove precursor peaks, 0: no, 1: yes, 2: yes + all charge reduced precursor peaks, 3: yes + precursor phosphate neutral loss peaks, default is '0'.", false, true),
    COMET_REMOVE_PRECURSOR_TOLERANCE("comet_remove_prec_tol", "Comet remove precursor tolerance, default is '1.5'.", false, true),
    COMET_CLEAR_MZ_RANGE_LOWER("comet_clear_mz_range_lower", "Comet clear mz range lower, default is '0.0'.", false, true),
    COMET_CLEAR_MZ_RANGE_UPPER("comet_clear_mz_range_upper", "Comet clear mz range upper, default is '0.0'.", false, true),
    COMET_ENZYME_TYPE("comet_enzyme_type", "Comet enzyme type, 1: semi-specific, 2: full-enzyme, 8: unspecific N-term, 9: unspecific C-term, default is '2'.", false, true),
    COMET_ISOTOPE_CORRECTION("comet_isotope_correction", "Comet isotope correction, 0: no correction, 1: 0, +1, 2: 0, +1, +2, 3: 0,+1,+2,+3, 4: -8,-4,0,+4,+8, default is '1'.", false, true),
    COMET_MIN_PREC_MASS("comet_min_prec_mass", "Comet minimum precursor mass, default is '600.0'.", false, true),
    COMET_MAX_PREC_MASS("comet_max_prec_mass", "Comet maximum precursor mass, default is '5000.0'.", false, true),
    COMET_MIN_PEP_LENGTH("comet_min_pep_length", "Comet minimum peptide length, default is '8'.", false, true),
    COMET_MAX_PEP_LENGTH("comet_max_pep_length", "Comet maximum peptide length, default is '30'.", false, true),
    COMET_MAX_FRAGMENT_CHARGE("comet_max_frag_charge", "Comet maximum fragment charge [1-5], default is '3'.", false, true),
    COMET_REMOVE_METH("comet_remove_meth", "Comet remove methionine, 1: true, 0: false, default is '0'.", false, true),
    COMET_BATCH_SIZE("comet_batch_size", "Comet spectrum batch size, '0' means load and search all spectra at once, default is '0'.", false, true),
    COMET_PTMS("comet_num_ptms", "Comet max number of variable PTMs per peptide, default is '10'.", false, true),
    COMET_REQ_PTMS("comet_req_ptms", "Comet require at least one variable PTM per peptide, 1: true, 0: false, default is '0'.", false, true),
    COMET_THEORETICAL_FRAGMENT_IONS("comet_theoretical_fragment_ions", "Comet theoretical_fragment_ions option, it is the correlation score type, 1: true, 0: false, default is '1'.", false, true),
    COMET_FRAGMENT_BIN_OFFSET("comet_frag_bin_offset", "Comet fragment bin offset, default is '0.01'.", false, true),
    COMET_NUM_MATCHES("comet_num_matches", "Comet maximum number of spectrum matches, default is '10'.", false, true),
    COMET_OUTPUT("comet_output", "Comet output type, PepXML, SQT, TXT or Percolator, default is 'PepXML'.", false, true),
    //////////////////////////////////
    // Tide specific parameters
    //////////////////////////////////
    TIDE_MIN_PEP_LENGTH("tide_min_pep_length", "Tide minumum peptide length, default is '6'.", false, true),
    TIDE_MAX_PEP_LENGTH("tide_max_pep_length", "Tide maximum peptide length, default is '30'.", false, true),
    TIDE_MIN_PREC_MASS("tide_min_prec_mass", "Tide minumum precursor mass, default is '200.0'.", false, true),
    TIDE_MAX_PREC_MASS("tide_max_prec_mass", "Tide maximum precursor mass, default is '7200.0'.", false, true),
    TIDE_MONOISOTOPIC("tide_monoisotopic", "Tide monoisotopic precursor, 1: true, 0: false, default is '1'.", false, true),
    TIDE_CLIP_N_TERM("tide_clip_n_term", "Tide clip n term methionine, 1: true, 0: false, default is '0'.", false, true),
    TIDE_PTMS("tide_num_ptms", "Tide max number of PTMs per peptide, default is no limit.", false, true),
    TIDE_PTMS_PER_TYPE("tide_num_ptms_per_type", "Tide max number of PTMs of each type per peptide, default is '2'.", false, true),
    TIDE_DIGESTION_TYPE("tide_digestion_type", "Tide digetion type (full-digest or partial-digest, true), default is 'full-digest'.", false, true),
    TIDE_PRINT_PEPTIDES("tide_print_peptides", "Tide print peptides, 1: true, 0: false, default is '0'.", false, true),
    TIDE_DECOY_FORMAT("tide_decoy_format", "Tide decoy fomat (none|shuffle|peptide-reverse|protein-reverse, true), default is 'none'.", false, true),
    TIDE_KEEP_TERM_AA("tide_keep_terminals", "Tide keep terminal amino acids when creating decoys (N|C|NC|none, true), default is 'NC'.", false, true),
    TIDE_DECOY_SEED("tide_decoy_seed", "Tide decoy seed, default is '1'.", false, true),
    TIDE_REMOVE_TEMP("tide_remove_temp", "Tide remove temp folders when the search is done, 1: true, 0: false, default is '1'.", false, true),
    TIDE_COMPUTE_P("tide_compute_p", "Tide compute exact p-values, 1: true, 0: false, default is '0'.", false, true),
    TIDE_COMPUTE_SP("tide_compute_sp", "Tide compute sp score, 1: true, 0: false, default is '0'.", false, true),
    TIDE_MIN_SPECTRUM_MZ("tide_min_spectrum_mz", "Tide minimum spectrum mz, default is '0.0'.", false, true),
    TIDE_MAX_SPECTRUM_MZ("tide_max_spectrum_mz", "Tide maximum spectrum mz, default is no limit.", false, true),
    TIDE_MIN_SPECTRUM_PEAKS("tide_min_spectrum_peaks", "Tide min spectrum peaks, default is '20'.", false, true),
    TIDE_SPECTRUM_CHARGES("tide_spectrum_charges", "Tide spectrum charges (1|2|3|all, true), default is 'all'.", false, true),
    TIDE_REMOVE_PREC("tide_remove_prec", "Tide remove precursor, 1: true, 0: false, default is '0'.", false, true),
    TIDE_REMOVE_PREC_TOL("tide_remove_prec_tol", "Tide remove precursor tolerance, default is '1.5'.", false, true),
    TIDE_USE_FLANKING("tide_use_flanking", "Tide use flanking peaks, 1: true, 0: false, default is '0'.", false, true),
    TIDE_USE_NEUTRAL_LOSSES("tide_use_neutral_losses", "Tide use neutral losses peaks, 1: true, 0: false, default is '0'.", false, true),
    TIDE_MZ_BIN_WIDTH("tide_mz_bin_width", "Tide mz bin width, default is '0.02'.", false, true),
    TIDE_MZ_BIN_OFFSET("tide_mz_bin_offset", "Tide mz bin offset, default is '0.0'.", false, true),
    TIDE_MAX_PSMS("tide_max_psms", "Tide maximum number of spectrum matches spectrum, default is '10'.", false, true),
    TIDE_EXPORT_TEXT("tide_export_text", "Tide export text file, 1: true, 0: false, default is '1'.", false, true),
    TIDE_EXPORT_SQT("tide_export_sqt", "Tide export SQT file, 1: true, 0: false, default is '0'.", false, true),
    TIDE_EXPORT_PEPXML("tide_export_pepxml", "Tide export pepxml, 1: true, 0: false, default is '0'.", false, true),
    TIDE_EXPORT_MZID("tide_export_mzid", "Tide export mzid, 1: true, 0: false, default is '0'.", false, true),
    TIDE_EXPORT_PIN("tide_export_pin", "Tide export Percolator input file, 1: true, 0: false, default is '0'.", false, true),
    TIDE_OUTPUT_FOLDER("tide_output_folder", "Tide output folder (relative to the Tide working folder, true), default is 'crux-output'.", false, true),
    TIDE_VERBOSITY("tide_verbosity", "Tide progress display verbosity (0|10|20|30|40|50|60, true), default is '30'.", false, true),
    TIDE_PROGRESS_INDICATOR("tide_progress_indicator", "Tide progress indicator frequency, default is '1000'.", false, true),
    TIDE_CONCAT("tide_concat", "Tide concatenate target and decoy results, 1: true, 0: false, default is '0'.", false, true),
    TIDE_STORE_SPECTRA("tide_store_spectra", "Tide file name in with to store the binary spectra, default is null, i.e., not set.", false, true),
    //////////////////////////////////
    // Andromeda specific parameters
    //////////////////////////////////
    ANDROMEDA_MAX_PEPTIDE_MASS("andromeda_max_pep_mass", "Andromeda maximum peptide mass, default is '4600.0'.", false, true),
    ANDROMEDA_MAX_COMBINATIONS("andromeda_max_comb", "Andromeda maximum combinations, default is '250'.", false, true),
    ANDROMEDA_TOP_PEAKS("andromeda_top_peaks", "Andromeda number of top peaks, default is '8'.", false, true),
    ANDROMEDA_TOP_PEAKS_WINDOW("andromeda_top_peaks_window", "Andromeda top peaks window width, default is '100'.", false, true),
    ANDROMEDA_INCL_WATER("andromeda_incl_water", "Andromeda account for water losses, 1: true, 0: false, default is '1'.", false, true),
    ANDROMEDA_INCL_AMMONIA("andromeda_incl_ammonia", "Andromeda account for ammonina losses, 1: true, 0: false, default is '1'.", false, true),
    ANDROMEDA_NEUTRAL_LOSSES("andromeda_neutral_losses", "Andromeda neutral losses are sequence dependent, 1: true, 0: false, default is '1'.", false, true),
    ANDROMEDA_FRAGMENT_ALL("andromeda_fragment_all", "Andromeda fragment all option, 1: true, 0: false, default is '0'.", false, true),
    ANDROMEDA_EMP_CORRECTION("andromeda_emp_correction", "Andromeda emperical correction, 1: true, 0: false, default is '1'.", false, true),
    ANDROMEDA_HIGHER_CHARGE("andromeda_higher_charge", "Andromeda higher charge option, 1: true, 0: false, default is '1'.", false, true),
    ANDROMEDA_FRAG_METHOD("andromeda_frag_method", "Andromeda fragmentation method, HCD, CID or EDT, default is 'CID'.", false, true),
    ANDROMEDA_MAX_MODS("andromeda_max_mods", "Andromeda maximum number of modifications, default is '5'.", false, true),
    ANDROMEDA_MIN_PEP_LENGTH("andromeda_min_pep_length", "Andromeda minimum peptide length when using no enzyme, default is '8'.", false, true),
    ANDROMEDA_MAX_PEP_LENGTH("andromeda_max_pep_length", "Andromeda maximum peptide length when using no enzyme, default is '30'.", false, true),
    ANDROMEDA_EQUAL_IL("andromeda_equal_il", "Andromeda whether I and L should be considered indistinguishable, 1: true, 0: false, default is '0'.", false, true),
    ANDROMEDA_MAX_PSMS("andromeda_max_psms", "Andromeda maximum number of spectrum matches spectrum, default is '10'.", false, true),
    ANDROMEDA_DECOY_MODE("andromeda_decoy_mode", "Andromeda decoy mode, none or decoy, default is 'none'.", false, true),
    //////////////////////////////////
    // MetaMorpheus specific parameters
    //////////////////////////////////
    META_MORPHEUS_MIN_PEP_LENGTH("meta_morpheus_min_pep_length", "MetaMorpheus minimum peptide length, default is '8'.", false, true),
    META_MORPHEUS_MAX_PEP_LENGTH("meta_morpheus_max_pep_length", "MetaMorpheus maximum peptide length, default is '30'.", false, true),
    META_MORPHEUS_SEARCH_TYPE("meta_morpheus_search_type", "MetaMorpheus search type, Classic, Modern or NonSpecific, default is 'Classic'.", false, true),
    META_MORPHEUS_NUM_PARTITIONS("meta_morpheus_num_partitions", "MetaMorpheus number of partitions when indexing, default is '1'.", false, true),
    META_MORPHEUS_DISSOCIATION_TYPE("meta_morpheus_dissociation_type", "MetaMorpheus dissociation type, HCD, CID, ECD or ETD, default is 'HCD'.", false, true),
    META_MORPHEUS_MAX_MODS_FOR_PEPTIDE("meta_morpheus_max_mods_for_peptide", "MetaMorpheus maximum modifications per peptide, default is '2'.", false, true),
    META_MORPHEUS_INITIATOR_METHIONINE_BEHAVIOR("meta_morpheus_meth", "MetaMorpheus initiator methionine behavior, Undefined, Retain, Cleave or Variable, default is 'Variable'.", false, true),
    META_MORPHEUS_SCORE_CUTOFF("meta_morpheus_score_cutoff", "MetaMorpheus score cutoff, default is '5.0'.", false, true),
    META_MORPHEUS_USE_DELTA_SCORE("meta_morpheus_use_delta_score", "MetaMorpheus use delta score, 1: true, 0: false, default is '0'.", false, true),
    META_MORPHEUS_FRAGMENTATION_TERMINUS("meta_morpheus_frag_term", "MetaMorpheus fragmentation terminus, Both, N or C, default is 'Both'.", false, true),
    META_MORPHEUS_MAX_FRAGMENTATION_SIZE("meta_morpheus_max_frag_size", "MetaMorpheus maximum fragmentation size, default is '30000'.", false, true),
    META_MORPHEUS_MASS_DIFF_ACCEPTOR_TYPE("meta_morpheus_mass_diff_acceptor_type", "MetaMorpheus mass difference acceptor type, Exact, OneMM, TwoMM, ThreeMM, PlusOrMinusThreeMM, ModOpen or Open, default is 'OneMM'.", false, true),
    META_MORPHEUS_WRITE_MZID("meta_morpheus_write_mzid", "MetaMorpheus write mzid, 1: true, 0: false, default is '1'.", false, true),
    META_MORPHEUS_WRITE_PEPXML("meta_morpheus_write_pepxml", "MetaMorpheus write pepxml, 1: true, 0: false, default is '0'.", false, true),
    META_MORPHEUS_USE_PROVIDED_PRECURSOR("meta_morpheus_use_provided_prec", "MetaMorpheus use provided precursor info, 1: true, 0: false, default is '1'.", false, true),
    META_MORPHEUS_DO_PREC_DECONVOLUTION("meta_morpheus_do_prec_deconv", "MetaMorpheus do precursor deconvolution, 1: true, 0: false, default is '1'.", false, true),
    META_MORPHEUS_DECONVOLUTION_INT_RATIO("meta_morpheus_deconv_int_ratio", "MetaMorpheus deconvolution intensity ratio, default is '3.0'.", false, true),
    META_MORPHEUS_DECONVOLUTION_MASS_TOL("meta_morpheus_deconv_mass_tol", "MetaMorpheus deoconvolution mass tolerance, default is '4.0'.", false, true),
    META_MORPHEUS_DECONVOLUTION_MASS_TOL_TYPE("meta_morpheus_deconv_mass_tol_type", "MetaMorpheus deoconvolution mass tolerance type, PPM or Absolute, default is 'PPM'.", false, true),
    META_MORPHEUS_TRIM_MS1_PEAKS("meta_morpheus_trim_ms1", "MetaMorpheus trim MS1 peaks, 1: true, 0: false, default is '0'.", false, true),
    META_MORPHEUS_TRIM_MSMS_PEAKS("meta_morpheus_trim_msms", "MetaMorpheus trim MSMS peaks, 1: true, 0: false, default is '1'.", false, true),
    META_MORPHEUS_NUM_PEAKS_PER_WINDOWS("meta_morpheus_num_peaks_per_window", "MetaMorpheus number of peaks per window, default is '200'.", false, true),
    META_MORPHEUS_MIN_ALLOWED_INT_RATIO_TO_BASE_PEAK("meta_morpheus_min_allowed_int_ratio_to_base_peak", "MetaMorpheus minium allowed intensity ratio to base peak, default is '0.01'.", false, true),
    META_MORPHEUS_WINDOW_WITH_THOMPSON("meta_morpheus_window_with_thompson", "MetaMorpheus window width in Thompson.", false, true),
    META_MORPHEUS_NUM_WINDOWS("meta_morpheus_num_windows", "MetaMorpheus number of windows", false, true),
    META_MORPHEUS_NORMALIZE_PEAKS_ACROSS_ALL_WINDOWS("meta_morpheus_norm_across_all_windows", "MetaMorpheus normalize peaks across all windows, 1: true, 0: false, default is '0'.", false, true),
    META_MORPHEUS_MOD_PEPTIDES_ARE_DIFFERENT("meta_morpheus_mod_peptides_are_different", "MetaMorpheus modified peptides are different, 1: true, 0: false, default is '0'.", false, true),
    META_MORPHEUS_NO_ONE_HIT_WONDERS("meta_morpheus_no_one_hit_wonders", "MetaMorpheus exclude one hit wonders, 1: true, 0: false, default is '0'.", false, true),
    META_MORPHEUS_SEARCH_TARGET("meta_morpheus_search_target", "MetaMorpheus search target sequences, 1: true, 0: false, default is '1'.", false, true),
    META_MORPHEUS_DECOY_TYPE("meta_morpheus_decoy_type", "MetaMorpheus decoy type, None, Reverse or Slide, default is 'None'.", false, true),
    META_MORPHEUS_MAX_MOD_ISOFORMS("meta_morpheus_max_mod_isoforms", "MetaMorpheus maximum modified isoforms, default is '1024'.", false, true),
    META_MORPHEUS_MIN_VARIANT_DEPTH("meta_morpheus_min_variant_depth", "MetaMorpheus minimum variant depth, default is '1'.", false, true),
    META_MORPHEUS_MAX_HETROZYGOUS_VARIANTS("meta_morpheus_max_hetrozygous_var", "MetaMorpheus maximum hetrozygous variants, default is '4'.", false, true),
    //////////////////////////////////
    // DirecTag specific parameters
    //////////////////////////////////
    DIRECTAG_TAG_LENGTH("directag_tag_length", "DirecTag tag length, default is '4'.", false, true),
    DIRECTAG_MAX_DYNAMIC_MODS("directag_max_var_mods", "DirecTag maximum variable modifications per sequence, default is '2'.", false, true),
    DIRECTAG_NUM_CHARGE_STATES("directag_charge_states", "DirecTag number of charge states considered, default is '3'.", false, true),
    DIRECTAG_DUPLICATE_SPECTRA("directag_duplicate_spectra", "DirecTag duplicate spectra per charge, 1: true, 0: false, default is '1'.", false, true),
    DIRECTAG_ISOTOPE_MZ_TOLERANCE("directag_isotope_tolerance", "DirecTag isotope mz tolerance, default is '0.25'.", false, true),
    DIRECTAG_DEISOTOPING_MODE("directag_deisotoping", "DirecTag deisotoping mode, default is '0', 0: no deisotoping, 1: precursor only, 2: precursor and candidate.", false, true),
    DIRECTAG_NUM_INTENSITY_CLASSES("directag_intensity_classes", "DirecTag number of intensity classses, default is '3'.", false, true),
    DIRECTAG_OUTPUT_SUFFIX("directag_output_suffix", "DirecTag output suffic, default is no suffix.", false, true),
    DIRECTAG_MAX_PEAK_COUNT("directag_max_peak_count", "DirecTag max peak count, default is '100'.", false, true),
    DIRECTAG_MAX_TAG_COUNT("directag_max_tag_count", "DirecTag maximum tag count, default is '10'.", false, true),
    DIRECTAG_TIC_CUTOFF_PERCENTAGE("directag_tic_cutoff", "DirecTag TIC cutoff in percent, default is '100'.", false, true),
    DIRECTAG_COMPLEMENT_MZ_TOLERANCE("directag_complement_tolerance", "DirecTag complement mz tolerance, default is '0.1'.", false, true),
    DIRECTAG_PRECURSOR_ADJUSTMENT_STEP("directag_adjustment_step", "DirecTag precursor adjustment step, default is '0.1'.", false, true),
    DIRECTAG_MIN_PRECURSOR_ADJUSTMENT("directag_min_adjustment", "DirecTag minimum precursor adjustment, default is '-0.5'.", false, true),
    DIRECTAG_MAX_PRECURSOR_ADJUSTMENT("directag_max_adjustment", "DirecTag maximum precursor adjustment, default is '1.5'.", false, true),
    DIRECTAG_INTENSITY_SCORE_WEIGHT("directag_intensity_weight", "DirecTag intensity score weight, default is '1.0'.", false, true),
    DIRECTAG_MZ_FIDELITY_SCORE_WEIGHT("directag_fidelity_weight", "DirecTag fidelity score weight, default is '1.0'.", false, true),
    DIRECTAG_COMPLEMENT_SCORE_WEIGHT("directag_complement_weight", "DirecTag complement_score_weight, default is '1.0'.", false, true),
    DIRECTAG_ADJUST_PRECURSOR_MASS("directag_adjust_precursor", "DirecTag adjust precursor, 1: true, 0: false, default is '0'.", false, true),
    DIRECTAG_USE_CHARGE_STATE_FROM_MS("directag_ms_charge_state", "DirecTag use charge state from M spectrum, 1: true, 0: false, default is '1'.", false, true),
    //////////////////////////////////
    // PepNovo+ specific parameters
    //////////////////////////////////
    PEPNOVO_HITLIST_LENGTH("pepnovo_hitlist_length", "PepNovo+ number of de novo solutions [0-2000], default is '10'.", false, true),
    PEPNOVO_ESTIMATE_CHARGE("pepnovo_estimate_charge", "PepNovo+ estimate precursor charge option. 1: true, 0: false, default is '1'.", false, true),
    PEPNOVO_CORRECT_PREC_MASS("pepnovo_correct_prec_mass", "PepNovo+ correct precursor mass option. 1: true, 0: false, default is '1'.", false, true),
    PEPNOVO_DISCARD_SPECTRA("pepnovo_discard_spectra", "PepNovo+ discard low quality spectra optoin. 1: true, 0: false, default is '1'.", false, true),
    PEPNOVO_FRAGMENTATION_MODEL("pepnovo_fragmentation_model", "PepNovo+ fragmentation model. Default is 'CID_IT_TRYP'.", false, true),
    PEPNOVO_GENERATE_BLAST("pepnovo_generate_blast", "PepNovo+ generate a BLAST query. 1: true, 0: false, default is '0'.", false, true),
    //////////////////////////////////
    // pNovo+ specific parameters
    //////////////////////////////////
    PNOVO_NUMBER_OF_PEPTIDES("pnovo_num_peptides", "pNovo+ number of peptides per spectrum, default is '10'.", false, true),
    PNOVO_LOWER_PRECURSOR_MASS("pnovo_lower_prec", "pNovo+ minimum precursor mass, default is '300'.", false, true),
    PNOVO_UPPER_PRECURSOR_MASS("pnovo_upper_prec", "pNovo+ maximum precursor mass, default is '5000'.", false, true),
    PNOVO_ACTIVATION_TYPE("pnovo_activation", "pNovo+ activation type (HCD, CID or EDT, true), default is 'HCD'.", false, true),
    //////////////////////////////////
    // Novor specific parameters
    //////////////////////////////////
    NOVOR_FRAGMENTATION("novor_fragmentation", "Novor fragmentation method, CID or HCD, default is 'HCD'.", false, true),
    NOVOR_MASS_ANALYZER("novor_mass_analyzer", "Novor mass analyzer, Trap, TOF, or FT, default is 'FT'.", false, true),
    //////////////////////////////////
    // Gene annotation preferences
    //////////////////////////////////
    USE_GENE_MAPPING("useGeneMapping", "If true gene mappings will be used and saved along with the project (UniProt databases only). (1: true, 0: false, default is '1')", false, true),
    UPDATE_GENE_MAPPING("updateGeneMapping", "If true gene mappings will be updated automatically from Ensembl (UniProt databases only). (1: true, 0: false, default is '1')", false, true),
    //////////////////////////////////
    // Spectrum annotation
    //////////////////////////////////
    ANNOTATION_LEVEL("annotation_level", "The intensity threshold to consider for annotation, e.g. using percentiles, 0.75 means that the 25% most intense peaks will be annotated, default is 0.75.", false, true),
    ANNOTATION_LEVEL_TYPE("annotation_level_type", "The type of the intensity threshold. (valid values: " + AnnotationParameters.IntensityThresholdType.getCommandLineOptions() + ", default is " + AnnotationParameters.IntensityThresholdType.percentile + ").", false, true),
    ANNOTATION_MZ_TOLERANCE("annotation_mz_tolerance", "The m/z tolerance to annotate peaks, default is equal to the search settings MS2 tolerance.", false, true),
    ANNOTATION_HIGH_RESOLUTION("annotation_high_resolution", "If true the most accurate peak will be selected within the m/z tolerance. (1: true, 0: false, default is '1')", false, true),
    //////////////////////////////////
    // Sequence matching
    //////////////////////////////////
    SEQUENCE_MATCHING_TYPE("sequence_matching_type", "The peptide to protein sequence matching type. (" + SequenceMatchingParameters.MatchingType.getCommandLineOptions()
            + ", default is " + SequenceMatchingParameters.MatchingType.indistiguishableAminoAcids + ")", false, true),
    SEQUENCE_MATCHING_X("sequence_matching_x", "The maximum share of X's in a sequence, 0.25 means 25% of X's, default is 0.25.", false, true),
    //////////////////////////////////
    // Import filters
    //////////////////////////////////
    IMPORT_PEPTIDE_LENGTH_MIN("import_peptide_length_min", "The minimal peptide length to consider when importing identification files, default is 8.", false, true),
    IMPORT_PEPTIDE_LENGTH_MAX("import_peptide_length_max", "The maximal peptide length to consider when importing identification files, default is 30.", false, true),
    IMPORT_MC_MIN("import_missed_cleavages_min", "The minimal number if missed cleavages to consider when importing identification files, default is no filter.", false, true),
    IMPORT_MC_MAX("import_missed_cleavages_max", "The maximal number if missed cleavages to consider when importing identification files, default is no filter.", false, true),
    IMPORT_PRECURSOR_MZ("import_precursor_mz", "The maximal precursor deviation to allow when importing identification files, the precursor tolerance by default.", false, true),
    IMPORT_PRECURSOR_MZ_PPM("import_precursor_mz_ppm", "Maximal precursor ion deviation unit: ppm (1) or Da (0), default is '1'.", false, true),
    EXCLUDE_UNKNOWN_PTMs("exclude_unknown_ptms", "If true peptides presenting unrecognized PTMs will be excluded. (1: true, 0: false, default is '1')", false, true),
    //////////////////////////////////
    // PTM localization parameters
    //////////////////////////////////
    PTM_SCORE("ptm_score", "The PTM probabilistic score to use for modification localization (" + ModificationLocalizationScore.getCommandLineOptions() + ", default is '1').", false, true),
    PTM_THRESHOLD("ptm_threshold", "The threshold to use for the modification localizatoin score. Default is 95.", false, true),
    SCORE_NEUTRAL_LOSSES("score_neutral_losses", "Include neutral losses in spectrum annotation of the PTM score (1: true, 0: false, default is '0').", false, true),
    PTM_SEQUENCE_MATCHING_TYPE("ptm_sequence_matching_type", "The modification to peptide sequence matching type. (" + SequenceMatchingParameters.MatchingType.getCommandLineOptions()
            + ", default is " + SequenceMatchingParameters.MatchingType.aminoAcid + ")", false, true),
    PTM_ALIGNMENT("ptm_alignment", "Align peptide ambiguously localized PTMs on confident sites (1: true, 0: false, default is '1').", false, true),
    //////////////////////////////////
    // Protein inference parameters
    //////////////////////////////////
    SIMPLIFY_GOUPS("simplify_groups", "Simplify protein groups, 1: yes, 0: no, default is 1.", false, true),
    SIMPLIFY_GOUPS_EVIDENCE("simplify_evidence", "Simplify protein groups based on the Uniprot protein evidence, 1: yes, 0: no, default is 1.", false, true),
    SIMPLIFY_GOUPS_CONFIDENCE("simplify_confidence", "Simplify protein groups based on the peptide confidence, 1: yes, 0: no, default is 1.", false, true),
    SIMPLIFY_GOUPS_CONFIDENCE_THRESHOLD("simplify_confidence_threshold", "Peptide confidence threshold below which a peptide is considered absent, default is 0.05.", false, true),
    SIMPLIFY_GOUPS_ENZYMATICITY("simplify_enzymaticity", "Simplify protein groups based on the peptide enzymaticity, 1: yes, 0: no, default is 1.", false, true),
    SIMPLIFY_GOUPS_VARIANT("simplify_variant", "Simplify protein groups based on the variant matching, 1: yes, 0: no, default is 1.", false, true),
    PROTEIN_INFERENCE_MODIFICATIONS("pi_modifications", "Account for modifications when mapping peptides to proteins, 1: yes, 0: no, default is 1.", false, true),
    //////////////////////////////////
    // Validation parameters
    //////////////////////////////////
    PSM_FDR("psm_fdr", "FDR at the PSM level in percent, default is 1.", false, true),
    PEPTIDE_FDR("peptide_fdr", "FDR at the peptide level in percent, default is 1.", false, true),
    PROTEIN_FDR("protein_fdr", "FDR at the protein level in percent, default is 1.", false, true),
    //////////////////////////////////
    // Fraction parameters
    //////////////////////////////////
    PROTEIN_FRACTION_MW_CONFIDENCE("protein_fraction_mw_confidence", "Minimum confidence required for a protein in the fraction MW plot (default 95%: '95.0').", false, true),
    //////////////////////////////////
    // FASTA parameters
    //////////////////////////////////
    FASTA_TARGET_DECOY("fasta_target_decoy", "FASTA file should be processed as target-decoy. 1: true, 0: false, default is '1'.", false, true),
    FASTA_DECOY_TAG("fasta_decoy_tag", "The flag for decoy proteins in the accession. Default is '-REVERSED'.", false, true),
    FASTA_DECOY_SUFFIX("fasta_decoy_type", "Decoy type. 0: prefix, 1: suffix, default is '1'.", false, true),
    FASTA_DECOY_FILE_TAG("fasta_decoy_file_tag", "The tag added after adding decoy sequences to a FASTA file. Default is '_concatenated_target_decoy'.", false, true);

    /**
     * Short Id for the CLI parameter.
     */
    public final String id;
    /**
     * Explanation for the CLI parameter.
     */
    public final String description;
    /**
     * Boolean indicating whether the parameter is mandatory.
     */
    public final boolean mandatory;
    /**
     * Boolean indicating whether this command line option needs an argument.
     */
    public final boolean hasArgument;

    /**
     * Private constructor managing the various variables for the enum
     * instances.
     *
     * @param id the id
     * @param description the description
     * @param mandatory is the parameter mandatory
     * @param hasArgument boolean indicating whether this command line option
     * needs an argument
     */
    private IdentificationParametersCLIParams(String id, String description, boolean mandatory, boolean hasArgument) {
        this.id = id;
        this.description = description;
        this.mandatory = mandatory;
        this.hasArgument = hasArgument;
    }

    /**
     * Creates the options for the command line interface based on the possible
     * values.
     *
     * @param aOptions the options object where the options will be added
     */
    public static void createOptionsCLI(Options aOptions) {
        for (IdentificationParametersCLIParams identificationParametersCLIParams : IdentificationParametersCLIParams.values()) {
            aOptions.addOption(identificationParametersCLIParams.id, identificationParametersCLIParams.hasArgument, identificationParametersCLIParams.description);
        }
    }

    /**
     * Returns the options as a string.
     *
     * @return the options as a string
     */
    public static String getOptionsAsString() {

        String output = "";

        output += "Parameters Files:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OUT.id) + " " + IdentificationParametersCLIParams.OUT.description + " (Mandatory)\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.IDENTIFICATION_PARAMETERS.id) + " " + IdentificationParametersCLIParams.IDENTIFICATION_PARAMETERS.description + " (Optional)\n";
        output += getParametersOptionsAsString();
        return output;
    }

    /**
     * Returns the options as a string.
     *
     * @return the options as a string
     */
    public static String getParametersOptionsAsString() {

        String output = "";

        output += "\n\nSearch Parameters:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PREC_TOL.id) + " " + IdentificationParametersCLIParams.PREC_TOL.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PREC_PPM.id) + " " + IdentificationParametersCLIParams.PREC_PPM.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.FRAG_TOL.id) + " " + IdentificationParametersCLIParams.FRAG_TOL.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIGESTION.id) + " " + IdentificationParametersCLIParams.DIGESTION.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ENZYME.id) + " " + IdentificationParametersCLIParams.ENZYME.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.SPECIFICITY.id) + " " + IdentificationParametersCLIParams.SPECIFICITY.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.FIXED_MODS.id) + " " + IdentificationParametersCLIParams.FIXED_MODS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.VARIABLE_MODS.id) + " " + IdentificationParametersCLIParams.VARIABLE_MODS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MIN_CHARGE.id) + " " + IdentificationParametersCLIParams.MIN_CHARGE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MAX_CHARGE.id) + " " + IdentificationParametersCLIParams.MAX_CHARGE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MC.id) + " " + IdentificationParametersCLIParams.MC.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.FI.id) + " " + IdentificationParametersCLIParams.FI.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.RI.id) + " " + IdentificationParametersCLIParams.RI.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MIN_ISOTOPE.id) + " " + IdentificationParametersCLIParams.MIN_ISOTOPE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MAX_ISOTOPE.id) + " " + IdentificationParametersCLIParams.MAX_ISOTOPE.description + "\n";

        output += "\n\nX!Tandem advanced parameters:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_DYNAMIC_RANGE.id) + " " + IdentificationParametersCLIParams.XTANDEM_DYNAMIC_RANGE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_NPEAKS.id) + " " + IdentificationParametersCLIParams.XTANDEM_NPEAKS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_MIN_FRAG_MZ.id) + " " + IdentificationParametersCLIParams.XTANDEM_MIN_FRAG_MZ.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_MIN_PEAKS.id) + " " + IdentificationParametersCLIParams.XTANDEM_MIN_PEAKS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_NOISE_SUPPRESSION.id) + " " + IdentificationParametersCLIParams.XTANDEM_NOISE_SUPPRESSION.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_MIN_PREC_MASS.id) + " " + IdentificationParametersCLIParams.XTANDEM_MIN_PREC_MASS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_PARENT_MONOISOTOPIC_MASS_ISOTOPE_ERROR.id) + " " + IdentificationParametersCLIParams.XTANDEM_PARENT_MONOISOTOPIC_MASS_ISOTOPE_ERROR.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_QUICK_ACETYL.id) + " " + IdentificationParametersCLIParams.XTANDEM_QUICK_ACETYL.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_QUICK_PYRO.id) + " " + IdentificationParametersCLIParams.XTANDEM_QUICK_PYRO.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_STP_BIAS.id) + " " + IdentificationParametersCLIParams.XTANDEM_STP_BIAS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_REFINE.id) + " " + IdentificationParametersCLIParams.XTANDEM_REFINE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_REFINE_EVALUE.id) + " " + IdentificationParametersCLIParams.XTANDEM_REFINE_EVALUE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_REFINE_UNANTICIPATED_CLEAVAGE.id) + " " + IdentificationParametersCLIParams.XTANDEM_REFINE_UNANTICIPATED_CLEAVAGE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_REFINE_SEMI.id) + " " + IdentificationParametersCLIParams.XTANDEM_REFINE_SEMI.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_REFINE_POTENTIAL_MOD_FULL_REFINEMENT.id) + " " + IdentificationParametersCLIParams.XTANDEM_REFINE_POTENTIAL_MOD_FULL_REFINEMENT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_REFINE_POINT_MUTATIONS.id) + " " + IdentificationParametersCLIParams.XTANDEM_REFINE_POINT_MUTATIONS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_REFINE_SNAPS.id) + " " + IdentificationParametersCLIParams.XTANDEM_REFINE_SNAPS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_REFINE_SPECTRUM_SYNTHESIS.id) + " " + IdentificationParametersCLIParams.XTANDEM_REFINE_SPECTRUM_SYNTHESIS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_EVALUE.id) + " " + IdentificationParametersCLIParams.XTANDEM_EVALUE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_OUTPUT_RESULTS.id) + " " + IdentificationParametersCLIParams.XTANDEM_OUTPUT_RESULTS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_OUTPUT_PROTEINS.id) + " " + IdentificationParametersCLIParams.XTANDEM_OUTPUT_PROTEINS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_OUTPUT_SEQUENCES.id) + " " + IdentificationParametersCLIParams.XTANDEM_OUTPUT_SEQUENCES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_OUTPUT_SPECTRA.id) + " " + IdentificationParametersCLIParams.XTANDEM_OUTPUT_SPECTRA.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_OUTPUT_HISTOGRAMS.id) + " " + IdentificationParametersCLIParams.XTANDEM_OUTPUT_HISTOGRAMS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.XTANDEM_SKYLINE.id) + " " + IdentificationParametersCLIParams.XTANDEM_SKYLINE.description + "\n";

        output += "\n\nMyriMatch advanced parameters:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_MIN_PEP_LENGTH.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_MIN_PEP_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_MAX_PEP_LENGTH.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_MAX_PEP_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_MIN_PREC_MASS.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_MIN_PREC_MASS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_MAX_PREC_MASS.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_MAX_PREC_MASS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_NUM_MATCHES.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_NUM_MATCHES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_PTMS.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_PTMS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_FRAGMENTATION.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_FRAGMENTATION.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_TERMINI.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_TERMINI.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_SMART_PLUS_THREE.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_SMART_PLUS_THREE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_XCORR.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_XCORR.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_TIC_CUTOFF.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_TIC_CUTOFF.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_INTENSTITY_CLASSES.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_INTENSTITY_CLASSES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_CLASS_MULTIPLIER.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_CLASS_MULTIPLIER.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_NUM_BATCHES.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_NUM_BATCHES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_MAX_PEAK_COUNT.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_MAX_PEAK_COUNT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MYRIMATCH_OUTPUT_FORMAT.id) + " " + IdentificationParametersCLIParams.MYRIMATCH_OUTPUT_FORMAT.description + "\n";

        output += "\n\nMS Amanda advanced parameters:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_DECOY.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_DECOY.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_DECOY_RANKING.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_DECOY_RANKING.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_INSTRUMENT.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_INSTRUMENT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_MAX_RANK.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_MAX_RANK.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_MONOISOTOPIC.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_MONOISOTOPIC.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_PERFORM_DEISOTOPING.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_PERFORM_DEISOTOPING.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_MAX_MOD.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_MAX_MOD.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_MAX_VAR_MOD.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_MAX_VAR_MOD.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_MAX_MOD_SITES.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_MAX_MOD_SITES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_MAX_NEUTRAL_LOSSES.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_MAX_NEUTRAL_LOSSES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_MAX_NEUTRAL_LOSSES_MODIFICATIONS.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_MAX_NEUTRAL_LOSSES_MODIFICATIONS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_MIN_PEPTIDE_LENGTH.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_MIN_PEPTIDE_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_MAX_PEPTIDE_LENGTH.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_MAX_PEPTIDE_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_LOADED_PROTEINS.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_LOADED_PROTEINS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_LOADED_SPECTRA.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_LOADED_SPECTRA.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MS_AMANDA_OUTPUT_FORMAT.id) + " " + IdentificationParametersCLIParams.MS_AMANDA_OUTPUT_FORMAT.description + "\n";

        output += "\n\nMS-GF+ advanced parameters:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MSGF_DECOY.id) + " " + IdentificationParametersCLIParams.MSGF_DECOY.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MSGF_INSTRUMENT.id) + " " + IdentificationParametersCLIParams.MSGF_INSTRUMENT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MSGF_FRAGMENTATION.id) + " " + IdentificationParametersCLIParams.MSGF_FRAGMENTATION.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MSGF_PROTOCOL.id) + " " + IdentificationParametersCLIParams.MSGF_PROTOCOL.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MSGF_TERMINI.id) + " " + IdentificationParametersCLIParams.MSGF_TERMINI.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MSGF_MIN_PEP_LENGTH.id) + " " + IdentificationParametersCLIParams.MSGF_MIN_PEP_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MSGF_MAX_PEP_LENGTH.id) + " " + IdentificationParametersCLIParams.MSGF_MAX_PEP_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MSGF_PTMS.id) + " " + IdentificationParametersCLIParams.MSGF_PTMS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MSGF_NUM_MATCHES.id) + " " + IdentificationParametersCLIParams.MSGF_NUM_MATCHES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MSGF_ADDITIONAL.id) + " " + IdentificationParametersCLIParams.MSGF_ADDITIONAL.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MSGF_TASKS.id) + " " + IdentificationParametersCLIParams.MSGF_TASKS.description + "\n";

        output += "\n\nOMSSA advanced parameters:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_LOW_INTENSITY.id) + " " + IdentificationParametersCLIParams.OMSSA_LOW_INTENSITY.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_HIGH_INTENSITY.id) + " " + IdentificationParametersCLIParams.OMSSA_HIGH_INTENSITY.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_INTENSITY_INCREMENT.id) + " " + IdentificationParametersCLIParams.OMSSA_INTENSITY_INCREMENT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_MIN_PEAKS.id) + " " + IdentificationParametersCLIParams.OMSSA_MIN_PEAKS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_REMOVE_PREC.id) + " " + IdentificationParametersCLIParams.OMSSA_REMOVE_PREC.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_ESTIMATE_CHARGE.id) + " " + IdentificationParametersCLIParams.OMSSA_ESTIMATE_CHARGE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_PLUS_ONE.id) + " " + IdentificationParametersCLIParams.OMSSA_PLUS_ONE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_MAX_FRACTION.id) + " " + IdentificationParametersCLIParams.OMSSA_MAX_FRACTION.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_PREC_PER_SPECTRUM.id) + " " + IdentificationParametersCLIParams.OMSSA_PREC_PER_SPECTRUM.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_SCALE_PREC.id) + " " + IdentificationParametersCLIParams.OMSSA_SCALE_PREC.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_SEQUENCES_IN_MEMORY.id) + " " + IdentificationParametersCLIParams.OMSSA_SEQUENCES_IN_MEMORY.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_METHIONINE.id) + " " + IdentificationParametersCLIParams.OMSSA_METHIONINE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_NEUTRON.id) + " " + IdentificationParametersCLIParams.OMSSA_NEUTRON.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_SINGLE_WINDOW_WIDTH.id) + " " + IdentificationParametersCLIParams.OMSSA_SINGLE_WINDOW_WIDTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_DOUBLE_WINDOW_WIDTH.id) + " " + IdentificationParametersCLIParams.OMSSA_DOUBLE_WINDOW_WIDTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_SINGLE_WINDOW_PEAKS.id) + " " + IdentificationParametersCLIParams.OMSSA_SINGLE_WINDOW_PEAKS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_DOUBLE_WINDOW_PEAKS.id) + " " + IdentificationParametersCLIParams.OMSSA_DOUBLE_WINDOW_PEAKS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_MIN_ANNOTATED_INTENSE_PEAKS.id) + " " + IdentificationParametersCLIParams.OMSSA_MIN_ANNOTATED_INTENSE_PEAKS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_MIN_ANNOTATED_PEAKS.id) + " " + IdentificationParametersCLIParams.OMSSA_MIN_ANNOTATED_PEAKS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_MAX_LADDERS.id) + " " + IdentificationParametersCLIParams.OMSSA_MAX_LADDERS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_MAX_FRAG_CHARGE.id) + " " + IdentificationParametersCLIParams.OMSSA_MAX_FRAG_CHARGE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_POSITIVE_IONS.id) + " " + IdentificationParametersCLIParams.OMSSA_POSITIVE_IONS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_FORWARD_IONS.id) + " " + IdentificationParametersCLIParams.OMSSA_FORWARD_IONS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_REWIND_IONS.id) + " " + IdentificationParametersCLIParams.OMSSA_REWIND_IONS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_MAX_FRAG_SERIES.id) + " " + IdentificationParametersCLIParams.OMSSA_MAX_FRAG_SERIES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_CORRELATION_CORRECTION.id) + " " + IdentificationParametersCLIParams.OMSSA_CORRELATION_CORRECTION.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_CONSECUTIVE_ION_PROBABILITY.id) + " " + IdentificationParametersCLIParams.OMSSA_CONSECUTIVE_ION_PROBABILITY.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_HITLIST_LENGTH_CHARGE.id) + " " + IdentificationParametersCLIParams.OMSSA_HITLIST_LENGTH_CHARGE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_ITERATIVE_SEQUENCE_EVALUE.id) + " " + IdentificationParametersCLIParams.OMSSA_ITERATIVE_SEQUENCE_EVALUE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_ITERATIVE_SPECTRUM_EVALUE.id) + " " + IdentificationParametersCLIParams.OMSSA_ITERATIVE_SPECTRUM_EVALUE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_ITERATIVE_REPLACE_EVALUE.id) + " " + IdentificationParametersCLIParams.OMSSA_ITERATIVE_REPLACE_EVALUE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_MIN_PEP_LENGTH.id) + " " + IdentificationParametersCLIParams.OMSSA_MIN_PEP_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_MAX_PEP_LENGTH.id) + " " + IdentificationParametersCLIParams.OMSSA_MAX_PEP_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_MAX_EVALUE.id) + " " + IdentificationParametersCLIParams.OMSSA_MAX_EVALUE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_HITLIST_LENGTH.id) + " " + IdentificationParametersCLIParams.OMSSA_HITLIST_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.OMSSA_FORMAT.id) + " " + IdentificationParametersCLIParams.OMSSA_FORMAT.description + "\n";

        output += "\n\nComet advanced parameters:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_MIN_PEAKS.id) + " " + IdentificationParametersCLIParams.COMET_MIN_PEAKS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_MIN_PEAK_INTENSITY.id) + " " + IdentificationParametersCLIParams.COMET_MIN_PEAK_INTENSITY.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_REMOVE_PRECURSOR.id) + " " + IdentificationParametersCLIParams.COMET_REMOVE_PRECURSOR.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_REMOVE_PRECURSOR_TOLERANCE.id) + " " + IdentificationParametersCLIParams.COMET_REMOVE_PRECURSOR_TOLERANCE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_CLEAR_MZ_RANGE_LOWER.id) + " " + IdentificationParametersCLIParams.COMET_CLEAR_MZ_RANGE_LOWER.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_CLEAR_MZ_RANGE_UPPER.id) + " " + IdentificationParametersCLIParams.COMET_CLEAR_MZ_RANGE_UPPER.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_ENZYME_TYPE.id) + " " + IdentificationParametersCLIParams.COMET_ENZYME_TYPE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_ISOTOPE_CORRECTION.id) + " " + IdentificationParametersCLIParams.COMET_ISOTOPE_CORRECTION.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_MIN_PREC_MASS.id) + " " + IdentificationParametersCLIParams.COMET_MIN_PREC_MASS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_MAX_PREC_MASS.id) + " " + IdentificationParametersCLIParams.COMET_MAX_PREC_MASS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_MIN_PEP_LENGTH.id) + " " + IdentificationParametersCLIParams.COMET_MIN_PEP_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_MAX_PEP_LENGTH.id) + " " + IdentificationParametersCLIParams.COMET_MAX_PEP_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_MAX_FRAGMENT_CHARGE.id) + " " + IdentificationParametersCLIParams.COMET_MAX_FRAGMENT_CHARGE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_REMOVE_METH.id) + " " + IdentificationParametersCLIParams.COMET_REMOVE_METH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_BATCH_SIZE.id) + " " + IdentificationParametersCLIParams.COMET_BATCH_SIZE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_PTMS.id) + " " + IdentificationParametersCLIParams.COMET_PTMS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_REQ_PTMS.id) + " " + IdentificationParametersCLIParams.COMET_REQ_PTMS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_THEORETICAL_FRAGMENT_IONS.id) + " " + IdentificationParametersCLIParams.COMET_THEORETICAL_FRAGMENT_IONS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_FRAGMENT_BIN_OFFSET.id) + " " + IdentificationParametersCLIParams.COMET_FRAGMENT_BIN_OFFSET.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_NUM_MATCHES.id) + " " + IdentificationParametersCLIParams.COMET_NUM_MATCHES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.COMET_OUTPUT.id) + " " + IdentificationParametersCLIParams.COMET_OUTPUT.description + "\n";

        output += "\n\nTide advanced parameters:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_MIN_PEP_LENGTH.id) + " " + IdentificationParametersCLIParams.TIDE_MIN_PEP_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_MAX_PEP_LENGTH.id) + " " + IdentificationParametersCLIParams.TIDE_MAX_PEP_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_MIN_PREC_MASS.id) + " " + IdentificationParametersCLIParams.TIDE_MIN_PREC_MASS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_MAX_PREC_MASS.id) + " " + IdentificationParametersCLIParams.TIDE_MAX_PREC_MASS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_MONOISOTOPIC.id) + " " + IdentificationParametersCLIParams.TIDE_MONOISOTOPIC.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_CLIP_N_TERM.id) + " " + IdentificationParametersCLIParams.TIDE_CLIP_N_TERM.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_PTMS.id) + " " + IdentificationParametersCLIParams.TIDE_PTMS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_PTMS_PER_TYPE.id) + " " + IdentificationParametersCLIParams.TIDE_PTMS_PER_TYPE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_DIGESTION_TYPE.id) + " " + IdentificationParametersCLIParams.TIDE_DIGESTION_TYPE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_PRINT_PEPTIDES.id) + " " + IdentificationParametersCLIParams.TIDE_PRINT_PEPTIDES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_DECOY_FORMAT.id) + " " + IdentificationParametersCLIParams.TIDE_DECOY_FORMAT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_KEEP_TERM_AA.id) + " " + IdentificationParametersCLIParams.TIDE_KEEP_TERM_AA.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_DECOY_SEED.id) + " " + IdentificationParametersCLIParams.TIDE_DECOY_SEED.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_REMOVE_TEMP.id) + " " + IdentificationParametersCLIParams.TIDE_REMOVE_TEMP.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_COMPUTE_P.id) + " " + IdentificationParametersCLIParams.TIDE_COMPUTE_P.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_COMPUTE_SP.id) + " " + IdentificationParametersCLIParams.TIDE_COMPUTE_SP.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_MIN_SPECTRUM_MZ.id) + " " + IdentificationParametersCLIParams.TIDE_MIN_SPECTRUM_MZ.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_MAX_SPECTRUM_MZ.id) + " " + IdentificationParametersCLIParams.TIDE_MAX_SPECTRUM_MZ.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_MIN_SPECTRUM_PEAKS.id) + " " + IdentificationParametersCLIParams.TIDE_MIN_SPECTRUM_PEAKS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_SPECTRUM_CHARGES.id) + " " + IdentificationParametersCLIParams.TIDE_SPECTRUM_CHARGES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_REMOVE_PREC.id) + " " + IdentificationParametersCLIParams.TIDE_REMOVE_PREC.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_REMOVE_PREC_TOL.id) + " " + IdentificationParametersCLIParams.TIDE_REMOVE_PREC_TOL.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_USE_FLANKING.id) + " " + IdentificationParametersCLIParams.TIDE_USE_FLANKING.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_USE_NEUTRAL_LOSSES.id) + " " + IdentificationParametersCLIParams.TIDE_USE_NEUTRAL_LOSSES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_MZ_BIN_WIDTH.id) + " " + IdentificationParametersCLIParams.TIDE_MZ_BIN_WIDTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_MZ_BIN_OFFSET.id) + " " + IdentificationParametersCLIParams.TIDE_MZ_BIN_OFFSET.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_MAX_PSMS.id) + " " + IdentificationParametersCLIParams.TIDE_MAX_PSMS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_EXPORT_TEXT.id) + " " + IdentificationParametersCLIParams.TIDE_EXPORT_TEXT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_EXPORT_SQT.id) + " " + IdentificationParametersCLIParams.TIDE_EXPORT_SQT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_EXPORT_PEPXML.id) + " " + IdentificationParametersCLIParams.TIDE_EXPORT_PEPXML.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_EXPORT_MZID.id) + " " + IdentificationParametersCLIParams.TIDE_EXPORT_MZID.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_EXPORT_PIN.id) + " " + IdentificationParametersCLIParams.TIDE_EXPORT_PIN.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_OUTPUT_FOLDER.id) + " " + IdentificationParametersCLIParams.TIDE_OUTPUT_FOLDER.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_VERBOSITY.id) + " " + IdentificationParametersCLIParams.TIDE_VERBOSITY.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_PROGRESS_INDICATOR.id) + " " + IdentificationParametersCLIParams.TIDE_PROGRESS_INDICATOR.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_CONCAT.id) + " " + IdentificationParametersCLIParams.TIDE_CONCAT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.TIDE_STORE_SPECTRA.id) + " " + IdentificationParametersCLIParams.TIDE_STORE_SPECTRA.description + "\n";

        output += "\n\nAndromeda advanced parameters:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_MAX_PEPTIDE_MASS.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_MAX_PEPTIDE_MASS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_MAX_COMBINATIONS.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_MAX_COMBINATIONS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_TOP_PEAKS.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_TOP_PEAKS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_TOP_PEAKS_WINDOW.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_TOP_PEAKS_WINDOW.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_INCL_WATER.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_INCL_WATER.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_INCL_AMMONIA.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_INCL_AMMONIA.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_NEUTRAL_LOSSES.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_NEUTRAL_LOSSES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_FRAGMENT_ALL.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_FRAGMENT_ALL.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_EMP_CORRECTION.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_EMP_CORRECTION.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_HIGHER_CHARGE.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_HIGHER_CHARGE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_FRAG_METHOD.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_FRAG_METHOD.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_MAX_MODS.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_MAX_MODS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_MIN_PEP_LENGTH.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_MIN_PEP_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_MAX_PEP_LENGTH.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_MAX_PEP_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_EQUAL_IL.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_EQUAL_IL.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_MAX_PSMS.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_MAX_PSMS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANDROMEDA_DECOY_MODE.id) + " " + IdentificationParametersCLIParams.ANDROMEDA_DECOY_MODE.description + "\n";

        output += "\n\nMetaMorpheus advanced parameters:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_MIN_PEP_LENGTH.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_MIN_PEP_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_MAX_PEP_LENGTH.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_MAX_PEP_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_SEARCH_TYPE.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_SEARCH_TYPE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_NUM_PARTITIONS.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_NUM_PARTITIONS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_DISSOCIATION_TYPE.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_DISSOCIATION_TYPE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_MAX_MODS_FOR_PEPTIDE.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_MAX_MODS_FOR_PEPTIDE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_INITIATOR_METHIONINE_BEHAVIOR.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_INITIATOR_METHIONINE_BEHAVIOR.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_SCORE_CUTOFF.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_SCORE_CUTOFF.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_USE_DELTA_SCORE.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_USE_DELTA_SCORE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_FRAGMENTATION_TERMINUS.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_FRAGMENTATION_TERMINUS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_MAX_FRAGMENTATION_SIZE.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_MAX_FRAGMENTATION_SIZE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_MASS_DIFF_ACCEPTOR_TYPE.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_MASS_DIFF_ACCEPTOR_TYPE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_WRITE_MZID.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_WRITE_MZID.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_WRITE_PEPXML.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_WRITE_PEPXML.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_USE_PROVIDED_PRECURSOR.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_USE_PROVIDED_PRECURSOR.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_DO_PREC_DECONVOLUTION.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_DO_PREC_DECONVOLUTION.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_DECONVOLUTION_INT_RATIO.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_DECONVOLUTION_INT_RATIO.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_DECONVOLUTION_MASS_TOL.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_DECONVOLUTION_MASS_TOL.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_DECONVOLUTION_MASS_TOL_TYPE.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_DECONVOLUTION_MASS_TOL_TYPE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_TRIM_MS1_PEAKS.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_TRIM_MS1_PEAKS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_TRIM_MSMS_PEAKS.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_TRIM_MSMS_PEAKS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_NUM_PEAKS_PER_WINDOWS.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_NUM_PEAKS_PER_WINDOWS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_MIN_ALLOWED_INT_RATIO_TO_BASE_PEAK.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_MIN_ALLOWED_INT_RATIO_TO_BASE_PEAK.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_WINDOW_WITH_THOMPSON.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_WINDOW_WITH_THOMPSON.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_NUM_WINDOWS.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_NUM_WINDOWS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_NORMALIZE_PEAKS_ACROSS_ALL_WINDOWS.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_NORMALIZE_PEAKS_ACROSS_ALL_WINDOWS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_MOD_PEPTIDES_ARE_DIFFERENT.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_MOD_PEPTIDES_ARE_DIFFERENT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_NO_ONE_HIT_WONDERS.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_NO_ONE_HIT_WONDERS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_SEARCH_TARGET.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_SEARCH_TARGET.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_DECOY_TYPE.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_DECOY_TYPE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_MAX_MOD_ISOFORMS.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_MAX_MOD_ISOFORMS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_MIN_VARIANT_DEPTH.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_MIN_VARIANT_DEPTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.META_MORPHEUS_MAX_HETROZYGOUS_VARIANTS.id) + " " + IdentificationParametersCLIParams.META_MORPHEUS_MAX_HETROZYGOUS_VARIANTS.description + "\n";

        
        output += "\n\nNovor:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.NOVOR_FRAGMENTATION.id) + " " + IdentificationParametersCLIParams.NOVOR_FRAGMENTATION.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.NOVOR_MASS_ANALYZER.id) + " " + IdentificationParametersCLIParams.NOVOR_MASS_ANALYZER.description + "\n";

        output += "\n\nPNovo+:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PNOVO_ACTIVATION_TYPE.id) + " " + IdentificationParametersCLIParams.PNOVO_ACTIVATION_TYPE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PNOVO_LOWER_PRECURSOR_MASS.id) + " " + IdentificationParametersCLIParams.PNOVO_LOWER_PRECURSOR_MASS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PNOVO_UPPER_PRECURSOR_MASS.id) + " " + IdentificationParametersCLIParams.PNOVO_UPPER_PRECURSOR_MASS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PNOVO_NUMBER_OF_PEPTIDES.id) + " " + IdentificationParametersCLIParams.PNOVO_NUMBER_OF_PEPTIDES.description + "\n";

        output += "\n\nPepNovo+:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PEPNOVO_HITLIST_LENGTH.id) + " " + IdentificationParametersCLIParams.PEPNOVO_HITLIST_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PEPNOVO_ESTIMATE_CHARGE.id) + " " + IdentificationParametersCLIParams.PEPNOVO_ESTIMATE_CHARGE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PEPNOVO_CORRECT_PREC_MASS.id) + " " + IdentificationParametersCLIParams.PEPNOVO_CORRECT_PREC_MASS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PEPNOVO_DISCARD_SPECTRA.id) + " " + IdentificationParametersCLIParams.PEPNOVO_DISCARD_SPECTRA.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PEPNOVO_GENERATE_BLAST.id) + " " + IdentificationParametersCLIParams.PEPNOVO_GENERATE_BLAST.description + "\n";

        output += "\n\nDirectTag advanced parameters:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_TAG_LENGTH.id) + " " + IdentificationParametersCLIParams.DIRECTAG_TAG_LENGTH.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_MAX_DYNAMIC_MODS.id) + " " + IdentificationParametersCLIParams.DIRECTAG_MAX_DYNAMIC_MODS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_NUM_CHARGE_STATES.id) + " " + IdentificationParametersCLIParams.DIRECTAG_NUM_CHARGE_STATES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_DUPLICATE_SPECTRA.id) + " " + IdentificationParametersCLIParams.DIRECTAG_DUPLICATE_SPECTRA.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_ISOTOPE_MZ_TOLERANCE.id) + " " + IdentificationParametersCLIParams.DIRECTAG_ISOTOPE_MZ_TOLERANCE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_DEISOTOPING_MODE.id) + " " + IdentificationParametersCLIParams.DIRECTAG_DEISOTOPING_MODE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_NUM_INTENSITY_CLASSES.id) + " " + IdentificationParametersCLIParams.DIRECTAG_NUM_INTENSITY_CLASSES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_OUTPUT_SUFFIX.id) + " " + IdentificationParametersCLIParams.DIRECTAG_OUTPUT_SUFFIX.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_MAX_PEAK_COUNT.id) + " " + IdentificationParametersCLIParams.DIRECTAG_MAX_PEAK_COUNT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_MAX_TAG_COUNT.id) + " " + IdentificationParametersCLIParams.DIRECTAG_MAX_TAG_COUNT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_TIC_CUTOFF_PERCENTAGE.id) + " " + IdentificationParametersCLIParams.DIRECTAG_TIC_CUTOFF_PERCENTAGE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_COMPLEMENT_MZ_TOLERANCE.id) + " " + IdentificationParametersCLIParams.DIRECTAG_COMPLEMENT_MZ_TOLERANCE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_PRECURSOR_ADJUSTMENT_STEP.id) + " " + IdentificationParametersCLIParams.DIRECTAG_PRECURSOR_ADJUSTMENT_STEP.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_MIN_PRECURSOR_ADJUSTMENT.id) + " " + IdentificationParametersCLIParams.DIRECTAG_MIN_PRECURSOR_ADJUSTMENT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_MAX_PRECURSOR_ADJUSTMENT.id) + " " + IdentificationParametersCLIParams.DIRECTAG_MAX_PRECURSOR_ADJUSTMENT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_INTENSITY_SCORE_WEIGHT.id) + " " + IdentificationParametersCLIParams.DIRECTAG_INTENSITY_SCORE_WEIGHT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_MZ_FIDELITY_SCORE_WEIGHT.id) + " " + IdentificationParametersCLIParams.DIRECTAG_MZ_FIDELITY_SCORE_WEIGHT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_COMPLEMENT_SCORE_WEIGHT.id) + " " + IdentificationParametersCLIParams.DIRECTAG_COMPLEMENT_SCORE_WEIGHT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_ADJUST_PRECURSOR_MASS.id) + " " + IdentificationParametersCLIParams.DIRECTAG_ADJUST_PRECURSOR_MASS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.DIRECTAG_USE_CHARGE_STATE_FROM_MS.id) + " " + IdentificationParametersCLIParams.DIRECTAG_USE_CHARGE_STATE_FROM_MS.description + "\n";

        output += "\n\nSpectrum Annotation:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANNOTATION_LEVEL_TYPE.id) + " " + IdentificationParametersCLIParams.ANNOTATION_LEVEL_TYPE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANNOTATION_LEVEL.id) + " " + IdentificationParametersCLIParams.ANNOTATION_LEVEL.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANNOTATION_MZ_TOLERANCE.id) + " " + IdentificationParametersCLIParams.ANNOTATION_MZ_TOLERANCE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.ANNOTATION_HIGH_RESOLUTION.id) + " " + IdentificationParametersCLIParams.ANNOTATION_HIGH_RESOLUTION.description + "\n";

        output += "\n\nSequence Matching:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.SEQUENCE_MATCHING_TYPE.id) + " " + IdentificationParametersCLIParams.SEQUENCE_MATCHING_TYPE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.SEQUENCE_MATCHING_X.id) + " " + IdentificationParametersCLIParams.SEQUENCE_MATCHING_X.description + "\n";

        output += "\n\nImport Filters:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.IMPORT_PEPTIDE_LENGTH_MIN.id) + " " + IdentificationParametersCLIParams.IMPORT_PEPTIDE_LENGTH_MIN.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.IMPORT_PEPTIDE_LENGTH_MAX.id) + " " + IdentificationParametersCLIParams.IMPORT_PEPTIDE_LENGTH_MAX.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.IMPORT_MC_MIN.id) + " " + IdentificationParametersCLIParams.IMPORT_MC_MIN.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.IMPORT_MC_MAX.id) + " " + IdentificationParametersCLIParams.IMPORT_MC_MAX.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.IMPORT_PRECURSOR_MZ.id) + " " + IdentificationParametersCLIParams.IMPORT_PRECURSOR_MZ.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.IMPORT_PRECURSOR_MZ_PPM.id) + " " + IdentificationParametersCLIParams.IMPORT_PRECURSOR_MZ_PPM.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.EXCLUDE_UNKNOWN_PTMs.id) + " " + IdentificationParametersCLIParams.EXCLUDE_UNKNOWN_PTMs.description + "\n";

        output += "\n\nPTM Localization:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PTM_SCORE.id) + " " + IdentificationParametersCLIParams.PTM_SCORE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PTM_THRESHOLD.id) + " " + IdentificationParametersCLIParams.PTM_THRESHOLD.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.SCORE_NEUTRAL_LOSSES.id) + " " + IdentificationParametersCLIParams.SCORE_NEUTRAL_LOSSES.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PTM_SEQUENCE_MATCHING_TYPE.id) + " " + IdentificationParametersCLIParams.PTM_SEQUENCE_MATCHING_TYPE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PTM_ALIGNMENT.id) + " " + IdentificationParametersCLIParams.PTM_ALIGNMENT.description + "\n";

        output += "\n\nGene Annotation:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.USE_GENE_MAPPING.id) + " " + IdentificationParametersCLIParams.USE_GENE_MAPPING.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.UPDATE_GENE_MAPPING.id) + " " + IdentificationParametersCLIParams.UPDATE_GENE_MAPPING.description + "\n";

        output += "\n\nProtein Inference:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.SIMPLIFY_GOUPS.id) + " " + IdentificationParametersCLIParams.SIMPLIFY_GOUPS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.SIMPLIFY_GOUPS_ENZYMATICITY.id) + " " + IdentificationParametersCLIParams.SIMPLIFY_GOUPS_ENZYMATICITY.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.SIMPLIFY_GOUPS_EVIDENCE.id) + " " + IdentificationParametersCLIParams.SIMPLIFY_GOUPS_EVIDENCE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.SIMPLIFY_GOUPS_CONFIDENCE.id) + " " + IdentificationParametersCLIParams.SIMPLIFY_GOUPS_CONFIDENCE.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.SIMPLIFY_GOUPS_CONFIDENCE_THRESHOLD.id) + " " + IdentificationParametersCLIParams.SIMPLIFY_GOUPS_CONFIDENCE_THRESHOLD.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.SIMPLIFY_GOUPS_VARIANT.id) + " " + IdentificationParametersCLIParams.SIMPLIFY_GOUPS_VARIANT.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PROTEIN_INFERENCE_MODIFICATIONS.id) + " " + IdentificationParametersCLIParams.PROTEIN_INFERENCE_MODIFICATIONS.description + "\n";

        output += "\n\nValidation Levels:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PSM_FDR.id) + " " + IdentificationParametersCLIParams.PSM_FDR.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PEPTIDE_FDR.id) + " " + IdentificationParametersCLIParams.PEPTIDE_FDR.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PROTEIN_FDR.id) + " " + IdentificationParametersCLIParams.PROTEIN_FDR.description + "\n";

        output += "\n\nFraction Analysis:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.PROTEIN_FRACTION_MW_CONFIDENCE.id) + " " + IdentificationParametersCLIParams.PROTEIN_FRACTION_MW_CONFIDENCE.description + "\n";

        output += "\n\nDatabase Processing:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.FASTA_TARGET_DECOY.id) + " " + IdentificationParametersCLIParams.FASTA_TARGET_DECOY.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.FASTA_DECOY_TAG.id) + " " + IdentificationParametersCLIParams.FASTA_DECOY_TAG.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.FASTA_DECOY_SUFFIX.id) + " " + IdentificationParametersCLIParams.FASTA_DECOY_SUFFIX.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.FASTA_DECOY_FILE_TAG.id) + " " + IdentificationParametersCLIParams.FASTA_DECOY_FILE_TAG.description + "\n";

        output += "\n\nHelp:\n\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.MODS.id) + " " + IdentificationParametersCLIParams.MODS.description + "\n";
        output += "-" + String.format(CommandLineUtils.formatter, IdentificationParametersCLIParams.USAGE.id) + " " + IdentificationParametersCLIParams.USAGE.description + "\n";

        return output;
    }
}
