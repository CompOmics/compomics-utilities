package com.compomics.util.test.experiment.io.identifications;

import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.identification.spectrum_assumptions.PeptideAssumption;
import com.compomics.util.experiment.io.identification.IdfileReader;
import com.compomics.util.experiment.io.identification.idfilereaders.InstaNovoIdfileReader;
import com.compomics.util.experiment.io.identification.idfilereaders.InstaNovoPlusIdfileReader;
import com.compomics.util.experiment.io.identification.idfilereaders.InstaNovoRefinedIdfileReader;
import com.compomics.util.experiment.mass_spectrometry.SpectrumProvider;
import com.compomics.util.experiment.mass_spectrometry.spectra.Precursor;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.parameters.identification.search.SearchParameters;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import junit.framework.TestCase;
import org.junit.Assert;

/**
 * Tests for InstaNovo v1.2.2 CSV readers.
 *
 * @author CompOmics
 */
public class TestInstaNovoIdfileReader extends TestCase {

    /**
     * Derived from the first row of the InstaNovo v1.2.2 transformer normalized
     * Zenodo sample file.
     */
    private static final String INSTANOVO_V1_2_2
            = "experiment_name,scan_number,spectrum_id,precursor_mz,precursor_charge,prediction_id,predictions,log_probs,token_log_probs,group,predictions_tokenised,delta_mass_ppm\n"
            + "SF_200217_U2OS_TiO2_HCD_OT_rep1,0,SF_200217_U2OS_TiO2_HCD_OT_rep1:0,419.314971923828,2,0,DM[UNIMOD:35]NS[UNIMOD:21]PK,-1147.98681640625,\"[-0.015801219269633293, -1.1395305395126343, -2.2013168334960938, -1.3749353885650635, -1.4705305099487305, -0.5675679445266724]\",no_group,\"D, M[UNIMOD:35], N, S[UNIMOD:21], P, K\",58846.475981092575\n";

    /**
     * Derived from the first row of the InstaNovo+ v1.2.2 no-refinement
     * normalized Zenodo sample file.
     */
    private static final String INSTANOVOPLUS_V1_2_2
            = "experiment_name,scan_number,spectrum_id,precursor_mz,precursor_charge,prediction_id,predictions,log_probs,token_log_probs,group,predictions_tokenised,delta_mass_ppm\n"
            + "SF_200217_U2OS_TiO2_HCD_OT_rep1,0,SF_200217_U2OS_TiO2_HCD_OT_rep1:0,419.314971923828,2,0,MC[UNIMOD:4]IPDQPM[UNIMOD:35]EVDNEDDAPLPPPEAR,-3.6934256553649902,,no_group,\"M, C[UNIMOD:4], I, P, D, Q, P, M[UNIMOD:35], E, V, D, N, E, D, D, A, P, L, P, P, P, E, A, R\",2282970.310323359\n";

    /**
     * Derived from the first row of the InstaNovo v1.2.2 combined refined
     * Zenodo sample file.
     */
    private static final String INSTANOVO_COMBINED_V1_2_2
            = "experiment_name,scan_number,spectrum_id,precursor_mz,precursor_charge,prediction_id,predictions,log_probs,token_log_probs,group,instanovo_predictions,instanovo_prediction_log_probability,instanovo_prediction_token_log_probabilities,instanovo_predictions_beam_0,instanovo_predictions_log_probability_beam_0,instanovo_predictions_token_log_probabilities_beam_0,instanovo_predictions_beam_1,instanovo_predictions_log_probability_beam_1,instanovo_predictions_token_log_probabilities_beam_1,instanovo_predictions_beam_2,instanovo_predictions_log_probability_beam_2,instanovo_predictions_token_log_probabilities_beam_2,instanovo_predictions_beam_3,instanovo_predictions_log_probability_beam_3,instanovo_predictions_token_log_probabilities_beam_3,instanovo_predictions_beam_4,instanovo_predictions_log_probability_beam_4,instanovo_predictions_token_log_probabilities_beam_4,instanovoplus_predictions,instanovoplus_prediction_log_probability,instanovoplus_prediction_token_log_probabilities,instanovoplus_unrefined_predictions,predictions_tokenised,delta_mass_ppm\n"
            + "SF_200217_U2OS_TiO2_HCD_OT_rep1,0,SF_200217_U2OS_TiO2_HCD_OT_rep1:0,419.314971923828,2,0,LIRPLLK,-0.6334811449050903,,no_group,\"['L', 'K', 'G', 'D', 'S[UNIMOD:21]', 'P', 'K']\",-10.102036476135254,\"[-1.716342806816101, -1.0499515533447266, -1.1343414783477783, -2.570066452026367, -1.3749353885650635, -1.4704134464263916, -0.5675679445266724]\",LKGDS[UNIMOD:21]PK,-10.102036476135254,\"[-1.716342806816101, -1.0499515533447266, -1.1343414783477783, -2.570066452026367, -1.3749353885650635, -1.4704134464263916, -0.5675679445266724]\",VKGDS[UNIMOD:21]PK,-11.082494735717773,\"[-2.8237648010253906, -1.0499515533447266, -1.1343414783477783, -2.570066452026367, -1.3749353885650635, -1.4704134464263916, -0.5675679445266724]\",SKGDS[UNIMOD:21]PK,-11.430251121520996,\"[-2.7461280822753906, -1.0499515533447266, -1.1343414783477783, -2.570066452026367, -1.3749353885650635, -1.4704134464263916, -0.5675679445266724]\",AKGDS[UNIMOD:21]PK,-11.492465019226074,\"[-3.1643409729003906, -1.0499515533447266, -1.1343414783477783, -2.570066452026367, -1.3749353885650635, -1.4704134464263916, -0.5675679445266724]\",PKGDS[UNIMOD:21]PK,-11.968438148498535,\"[-2.6694679260253906, -1.0499515533447266, -1.1343414783477783, -2.570066452026367, -1.3749353885650635, -1.4704134464263916, -0.5675679445266724]\",\"['L', 'I', 'R', 'P', 'L', 'L', 'K']\",-0.6334811449050903,,\"['L', 'K', 'G', 'D', 'S[UNIMOD:21]', 'P', 'K']\",\"L, I, R, P, L, L, K\",17862.82765389216\n";

    /**
     * Tests registration and parsing of the three supported InstaNovo modes.
     *
     * @throws Exception if an exception occurs
     */
    public void testInstaNovoReaders() throws Exception {

        Assert.assertNotNull(Advocate.getAdvocate("InstaNovo"));
        Assert.assertNotNull(Advocate.getAdvocate("InstaNovo+"));
        Assert.assertNotNull(Advocate.getAdvocate("InstaNovo with refinement"));
        SimpleSpectrumProvider spectrumProvider = new SimpleSpectrumProvider();
        SearchParameters searchParameters = new SearchParameters();

        assertReader("test.instanovo.csv", Advocate.instanovo.getIndex(), spectrumProvider, searchParameters);
        assertReader("test.instanovoplus.csv", Advocate.instanovoPlus.getIndex(), spectrumProvider, searchParameters);
        assertReader("test.instanovo.refined.csv", Advocate.instanovoRefined.getIndex(), spectrumProvider, searchParameters);
    }

    /**
     * Tests service registration for the three InstaNovo readers.
     *
     * @throws Exception if an exception occurs
     */
    public void testInstaNovoReaderServiceRegistration() throws Exception {

        InputStream serviceStream = getClass().getClassLoader().getResourceAsStream(
                "META-INF/services/com.compomics.util.experiment.io.identification.IdfileReader"
        );

        Assert.assertNotNull(serviceStream);

        byte[] bytes = new byte[serviceStream.available()];
        serviceStream.read(bytes);

        String serviceFile = new String(bytes, StandardCharsets.UTF_8);

        Assert.assertTrue(serviceFile.contains(InstaNovoIdfileReader.class.getName()));
        Assert.assertTrue(serviceFile.contains(InstaNovoPlusIdfileReader.class.getName()));
        Assert.assertTrue(serviceFile.contains(InstaNovoRefinedIdfileReader.class.getName()));
    }

    /**
     * Tests invalid headers.
     *
     * @throws Exception if an exception occurs
     */
    public void testMissingColumns() throws Exception {

        File csvFile = writeCsv("missing.instanovo.csv", "experiment_name,scan_number,predictions\nexample,0,PEPTIDE\n");
        IdfileReader idfileReader = new InstaNovoIdfileReader(csvFile);

        try {
            idfileReader.getAllSpectrumMatches(new SimpleSpectrumProvider(), null, new SearchParameters());
            Assert.fail("Expected invalid InstaNovo CSV columns to fail.");
        } catch (IllegalArgumentException e) {
            Assert.assertTrue(e.getMessage().contains("Mandatory"));
        }
    }

    /**
     * Tests parsing rows derived from the InstaNovo v1.2.2 sample files.
     *
     * @throws Exception if an exception occurs
     */
    public void testInstaNovoVersion122SampleRows() throws Exception {

        assertSampleReader(
                new InstaNovoIdfileReader(writeCsv("sample.instanovo.csv", INSTANOVO_V1_2_2)),
                Advocate.instanovo.getIndex(),
                "DMNSPK",
                2
        );

        assertSampleReader(
                new InstaNovoPlusIdfileReader(writeCsv("sample.instanovoplus.csv", INSTANOVOPLUS_V1_2_2)),
                Advocate.instanovoPlus.getIndex(),
                "MCIPDQPMEVDNEDDAPLPPPEAR",
                2
        );

        assertSampleReader(
                new InstaNovoRefinedIdfileReader(writeCsv("sample.instanovo.refined.csv", INSTANOVO_COMBINED_V1_2_2)),
                Advocate.instanovoRefined.getIndex(),
                "LIRPLLK",
                0
        );
    }

    /**
     * Tests matching realistic spectrum titles by scan tokens without positional
     * scan-number fallback.
     *
     * @throws Exception if an exception occurs
     */
    public void testSpectrumTitleLookupWithRealisticTitles() throws Exception {

        File csvFile = writeCsv(
                "realistic-titles.instanovo.csv",
                "experiment_name,scan_number,spectrum_id,precursor_mz,precursor_charge,prediction_id,predictions,log_probs\n"
                + "example,1,example:1,419.314971923828,2,0,PEPTIDE,-1.0\n"
        );

        IdfileReader idfileReader = new InstaNovoIdfileReader(csvFile);
        SimpleSpectrumProvider spectrumProvider = new SimpleSpectrumProvider(
                new String[]{"example"},
                new String[]{"controllerType=0 controllerNumber=1 scan=1", "controllerType=0 controllerNumber=1 scan=2"}
        );
        ArrayList<SpectrumMatch> spectrumMatches = idfileReader.getAllSpectrumMatches(spectrumProvider, null, new SearchParameters());

        Assert.assertEquals(1, spectrumMatches.size());
        Assert.assertEquals("controllerType=0 controllerNumber=1 scan=1", spectrumMatches.get(0).getSpectrumTitle());
    }

    /**
     * Tests matching InstaNovo positional spectrum ids to descriptive MGF
     * titles.
     *
     * @throws Exception if an exception occurs
     */
    public void testSpectrumTitleLookupWithPositionalSpectrumId() throws Exception {

        File csvFile = writeCsv(
                "positional-titles.instanovo.csv",
                "experiment_name,scan_number,spectrum_id,precursor_mz,precursor_charge,prediction_id,predictions,log_probs\n"
                + "example,0,example:0,419.314971923828,2,0,PEPTIDE,-1.0\n"
        );

        IdfileReader idfileReader = new InstaNovoIdfileReader(csvFile);
        SimpleSpectrumProvider spectrumProvider = new SimpleSpectrumProvider(
                new String[]{"example"},
                new String[]{"Cmpd 3543, +MSn(450.6095), 22.5 min", "Cmpd 3544, +MSn(697.8400), 22.5 min"}
        );
        ArrayList<SpectrumMatch> spectrumMatches = idfileReader.getAllSpectrumMatches(spectrumProvider, null, new SearchParameters());

        Assert.assertEquals(1, spectrumMatches.size());
        Assert.assertEquals("Cmpd 3543, +MSn(450.6095), 22.5 min", spectrumMatches.get(0).getSpectrumTitle());
    }

    /**
     * Tests charge parsing robustness.
     *
     * @throws Exception if an exception occurs
     */
    public void testChargeParsingSkipsInvalidRows() throws Exception {

        File csvFile = writeCsv(
                "charges.instanovo.csv",
                "experiment_name,scan_number,spectrum_id,precursor_mz,precursor_charge,prediction_id,predictions,log_probs\n"
                + "example,0,example:0,419.314971923828,not-a-charge,0,PEPTIDE,-1.0\n"
                + "example,1,example:1,419.314971923828, 2 ,0,PEPTIDE,-1.0\n"
        );

        IdfileReader idfileReader = new InstaNovoIdfileReader(csvFile);
        ArrayList<SpectrumMatch> spectrumMatches = idfileReader.getAllSpectrumMatches(new SimpleSpectrumProvider(), null, new SearchParameters());

        Assert.assertEquals(1, spectrumMatches.size());
        Assert.assertEquals("1", spectrumMatches.get(0).getSpectrumTitle());
    }

    /**
     * Tests all UniMod annotations from the InstaNovo v1.2.2 default residue
     * configuration.
     *
     * @throws Exception if an exception occurs
     */
    public void testDefaultInstaNovoModifications() throws Exception {

        String header = "experiment_name,scan_number,spectrum_id,precursor_mz,precursor_charge,prediction_id,predictions,log_probs,token_log_probs,group,predictions_tokenised,delta_mass_ppm\n";
        File csvFile = writeCsv(
                "default-modifications.instanovo.csv",
                header
                + "sample,0,sample:0,419.314971923828,2,0,M[UNIMOD:35]C[UNIMOD:4]N[UNIMOD:7]Q[UNIMOD:7]R[UNIMOD:7]P[UNIMOD:35]S[UNIMOD:21]T[UNIMOD:21]Y[UNIMOD:21],-1.0,,no_group,,0.0\n"
                + "sample,1,sample:1,419.314971923828,2,0,[UNIMOD:1]ACD,-1.0,,no_group,,0.0\n"
                + "sample,2,sample:2,419.314971923828,2,0,[UNIMOD:5]ACD,-1.0,,no_group,,0.0\n"
                + "sample,3,sample:3,419.314971923828,2,0,[UNIMOD:385]CPEP,-1.0,,no_group,,0.0\n"
                + "sample,4,sample:4,419.314971923828,2,0,[UNIMOD:385]NPEP,-1.0,,no_group,,0.0\n"
        );

        IdfileReader idfileReader = new InstaNovoIdfileReader(csvFile);
        ArrayList<SpectrumMatch> spectrumMatches = idfileReader.getAllSpectrumMatches(new SimpleSpectrumProvider(), null, new SearchParameters());

        Assert.assertEquals(5, spectrumMatches.size());

        PeptideAssumption residueModifiedAssumption = getFirstAssumption(spectrumMatches, "0", Advocate.instanovo.getIndex());

        Assert.assertEquals("MCNQRPSTY", residueModifiedAssumption.getPeptide().getSequence());
        assertModification(residueModifiedAssumption, "Oxidation of M", 1);
        assertModification(residueModifiedAssumption, "Carbamidomethylation of C", 2);
        assertModification(residueModifiedAssumption, "Deamidation of N", 3);
        assertModification(residueModifiedAssumption, "Deamidation of Q", 4);
        assertModification(residueModifiedAssumption, "Citrullination of R", 5);
        assertModification(residueModifiedAssumption, "Oxidation of P", 6);
        assertModification(residueModifiedAssumption, "Phosphorylation of S", 7);
        assertModification(residueModifiedAssumption, "Phosphorylation of T", 8);
        assertModification(residueModifiedAssumption, "Phosphorylation of Y", 9);

        assertModification(getFirstAssumption(spectrumMatches, "1", Advocate.instanovo.getIndex()), "Acetylation of peptide N-term", 0);
        assertModification(getFirstAssumption(spectrumMatches, "2", Advocate.instanovo.getIndex()), "Carbamilation of protein N-term", 0);
        assertModification(getFirstAssumption(spectrumMatches, "3", Advocate.instanovo.getIndex()), "Pyrolidone from carbamidomethylated C", 1);
        assertModification(getFirstAssumption(spectrumMatches, "4", Advocate.instanovo.getIndex()), "Ammonia loss from N", 1);
    }

    /**
     * Asserts one reader.
     *
     * @param fileName the file name
     * @param advocateIndex the expected advocate index
     * @param spectrumProvider the spectrum provider
     * @param searchParameters the search parameters
     *
     * @throws Exception if an exception occurs
     */
    private void assertReader(
            String fileName,
            int advocateIndex,
            SpectrumProvider spectrumProvider,
            SearchParameters searchParameters
    ) throws Exception {

        File csvFile = writeCsv(
                fileName,
                "experiment_name,scan_number,spectrum_id,precursor_mz,precursor_charge,prediction_id,predictions,log_probs,token_log_probs,group,predictions_tokenised,delta_mass_ppm\n"
                + "example,0,example:0,419.314971923828,2,0,DM[UNIMOD:35]NS[UNIMOD:21]PK,-10.0,\"[-1.0]\",no_group,\"D, M[UNIMOD:35], N, S[UNIMOD:21], P, K\",0.0\n"
        );

        IdfileReader idfileReader;
        if (fileName.endsWith(InstaNovoPlusIdfileReader.EXTENSION)) {
            idfileReader = new InstaNovoPlusIdfileReader(csvFile);
        } else if (fileName.endsWith(InstaNovoRefinedIdfileReader.EXTENSION)) {
            idfileReader = new InstaNovoRefinedIdfileReader(csvFile);
        } else {
            idfileReader = new InstaNovoIdfileReader(csvFile);
        }

        Assert.assertNotNull(idfileReader);

        ArrayList<SpectrumMatch> spectrumMatches = idfileReader.getAllSpectrumMatches(spectrumProvider, null, searchParameters);
        Assert.assertEquals(1, spectrumMatches.size());
        SpectrumMatch spectrumMatch = spectrumMatches.get(0);
        Assert.assertEquals("example", spectrumMatch.getSpectrumFile());
        Assert.assertEquals("0", spectrumMatch.getSpectrumTitle());

        TreeMap<Double, ArrayList<PeptideAssumption>> assumptions = spectrumMatch.getAllPeptideAssumptions(advocateIndex);
        Assert.assertNotNull(assumptions);
        PeptideAssumption peptideAssumption = assumptions.firstEntry().getValue().get(0);
        Assert.assertEquals("DMNSPK", peptideAssumption.getPeptide().getSequence());
        Assert.assertEquals(2, peptideAssumption.getPeptide().getVariableModifications().length);
        Assert.assertTrue(idfileReader.getSoftwareVersions().containsKey(Advocate.getAdvocate(advocateIndex).getName()));

        if (advocateIndex == Advocate.instanovoRefined.getIndex()) {
            Assert.assertTrue(idfileReader.getSoftwareVersions().containsKey(Advocate.instanovo.getName()));
            Assert.assertTrue(idfileReader.getSoftwareVersions().containsKey(Advocate.instanovoPlus.getName()));
        }
    }

    /**
     * Asserts a reader using sample v1.2.2 CSV content.
     *
     * @param idfileReader the reader
     * @param advocateIndex the expected advocate index
     * @param expectedSequence the expected peptide sequence
     * @param expectedVariableModifications the expected number of variable
     * modifications
     *
     * @throws Exception if an exception occurs
     */
    private void assertSampleReader(
            IdfileReader idfileReader,
            int advocateIndex,
            String expectedSequence,
            int expectedVariableModifications
    ) throws Exception {

        ArrayList<SpectrumMatch> spectrumMatches = idfileReader.getAllSpectrumMatches(new SimpleSpectrumProvider(), null, new SearchParameters());

        Assert.assertEquals(1, spectrumMatches.size());

        SpectrumMatch spectrumMatch = spectrumMatches.get(0);

        Assert.assertEquals("SF_200217_U2OS_TiO2_HCD_OT_rep1", spectrumMatch.getSpectrumFile());
        Assert.assertEquals("0", spectrumMatch.getSpectrumTitle());

        TreeMap<Double, ArrayList<PeptideAssumption>> assumptions = spectrumMatch.getAllPeptideAssumptions(advocateIndex);

        Assert.assertNotNull(assumptions);

        PeptideAssumption peptideAssumption = assumptions.firstEntry().getValue().get(0);

        Assert.assertEquals(expectedSequence, peptideAssumption.getPeptide().getSequence());
        Assert.assertEquals(expectedVariableModifications, peptideAssumption.getPeptide().getVariableModifications().length);
        Assert.assertTrue(idfileReader.getSoftwareVersions().containsKey(Advocate.getAdvocate(advocateIndex).getName()));

        if (advocateIndex == Advocate.instanovoRefined.getIndex()) {
            Assert.assertTrue(idfileReader.getSoftwareVersions().containsKey(Advocate.instanovo.getName()));
            Assert.assertTrue(idfileReader.getSoftwareVersions().containsKey(Advocate.instanovoPlus.getName()));
        }
    }

    /**
     * Returns the first assumption for a spectrum title.
     *
     * @param spectrumMatches the spectrum matches
     * @param spectrumTitle the spectrum title
     * @param advocateIndex the advocate index
     *
     * @return the first peptide assumption
     */
    private PeptideAssumption getFirstAssumption(ArrayList<SpectrumMatch> spectrumMatches, String spectrumTitle, int advocateIndex) {

        for (SpectrumMatch spectrumMatch : spectrumMatches) {

            if (spectrumMatch.getSpectrumTitle().equals(spectrumTitle)) {

                TreeMap<Double, ArrayList<PeptideAssumption>> assumptions = spectrumMatch.getAllPeptideAssumptions(advocateIndex);

                Assert.assertNotNull(assumptions);

                return assumptions.firstEntry().getValue().get(0);
            }
        }

        Assert.fail("No spectrum match found for title " + spectrumTitle + ".");

        return null;
    }

    /**
     * Asserts a modification match.
     *
     * @param peptideAssumption the peptide assumption
     * @param modification the modification name
     * @param site the modification site
     */
    private void assertModification(PeptideAssumption peptideAssumption, String modification, int site) {

        for (ModificationMatch modificationMatch : peptideAssumption.getPeptide().getVariableModifications()) {

            if (modificationMatch.getModification().equals(modification) && modificationMatch.getSite() == site) {
                return;
            }
        }

        Assert.fail("Modification " + modification + " at site " + site + " not found.");
    }

    /**
     * Writes a temporary CSV file.
     *
     * @param fileName the file name
     * @param content the content
     *
     * @return the CSV file
     *
     * @throws IOException if an IOException occurs
     */
    private File writeCsv(String fileName, String content) throws IOException {

        File file = File.createTempFile(fileName, "");
        file.deleteOnExit();

        try (FileWriter writer = new FileWriter(file)) {
            writer.write(content);
        }

        return file;
    }

    /**
     * Simple spectrum provider for tests.
     */
    private static class SimpleSpectrumProvider implements SpectrumProvider {

        /**
         * File names without extensions.
         */
        private final String[] fileNames;
        /**
         * Spectrum titles.
         */
        private final String[] titles;

        /**
         * Default constructor.
         */
        private SimpleSpectrumProvider() {
            this(new String[]{"example"}, new String[]{"0", "1", "2", "3", "4"});
        }

        /**
         * Constructor.
         *
         * @param fileNames the file names
         * @param titles the spectrum titles
         */
        private SimpleSpectrumProvider(String[] fileNames, String[] titles) {
            this.fileNames = fileNames;
            this.titles = titles;
        }

        @Override
        public Spectrum getSpectrum(String fileNameWithoutExtension, String spectrumTitle) {
            return null;
        }

        @Override
        public Precursor getPrecursor(String fileNameWithoutExtension, String spectrumTitle) {
            return null;
        }

        @Override
        public ArrayList<String> getPostcursorSpectrumTitles(String fileNameWithoutExtension, String spectrumTitle) {
            return null;
        }

        @Override
        public double getPrecursorMz(String fileNameWithoutExtension, String spectrumTitle) {
            return 0;
        }

        @Override
        public double getPrecursorRt(String fileNameWithoutExtension, String spectrumTitle) {
            return 0;
        }

        @Override
        public int getSpectrumLevel(String fileNameWithoutExtension, String spectrumTitle) {
            return 2;
        }

        @Override
        public double[][] getPeaks(String fileNameWithoutExtension, String spectrumTitle) {
            return null;
        }

        @Override
        public double getMinPrecMz(String fileNameWithoutExtension) {
            return 0;
        }

        @Override
        public double getMaxPrecMz(String fileNameWithoutExtension) {
            return 0;
        }

        @Override
        public double getMaxPrecInt(String fileNameWithoutExtension) {
            return 0;
        }

        @Override
        public double getMaxPrecRT(String fileNameWithoutExtension) {
            return 0;
        }

        @Override
        public double getMinPrecMz() {
            return 0;
        }

        @Override
        public double getMaxPrecMz() {
            return 0;
        }

        @Override
        public double getMaxPrecInt() {
            return 0;
        }

        @Override
        public double getMaxPrecRT() {
            return 0;
        }

        @Override
        public String[] getOrderedFileNamesWithoutExtensions() {
            return fileNames;
        }

        @Override
        public String[] getSpectrumTitles(String fileNameWithoutExtension) {
            return titles;
        }

        @Override
        public HashMap<String, String> getFilePaths() {
            return new HashMap<>();
        }

        @Override
        public HashMap<String, String> getCmsFilePaths() {
            return new HashMap<>();
        }

        @Override
        public void close() {
            // Nothing to close.
        }
    }
}
