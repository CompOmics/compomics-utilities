package com.compomics.util.experiment.biology.genes;

import com.compomics.util.Util;
import com.compomics.util.experiment.biology.genes.ensembl.EnsemblVersion;
import com.compomics.util.experiment.biology.genes.ensembl.GeneMapping;
import com.compomics.util.experiment.biology.genes.go.GoMapping;
import com.compomics.util.experiment.biology.taxonomy.SpeciesFactory;
import com.compomics.util.experiment.biology.taxonomy.mappings.EnsemblGenomesSpecies.EnsemblGenomeDivision;
import com.compomics.util.experiment.io.biology.protein.FastaSummary;
import com.compomics.util.gui.waiting.waitinghandlers.ProgressDialogX;
import com.compomics.util.parameters.identification.advanced.GeneParameters;
import com.compomics.util.experiment.io.biology.protein.ProteinDetailsProvider;
import com.compomics.util.experiment.io.biology.protein.SequenceProvider;
import com.compomics.util.io.IoUtil;
import com.compomics.util.waiting.WaitingHandler;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 * Class used to map proteins to gene information.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class ProteinGeneDetailsProvider {

    /**
     * The separator used to separate line contents.
     */
    public final static String SEPARATOR = "\t";
    /**
     * The folder where gene mapping files are stored.
     */
    private static String GENE_MAPPING_FOLDER = System.getProperty("user.home") + "/.compomics/gene_mappings/";
    /**
     * The subfolder relative to the jar file where gene mapping files are
     * stored in tools.
     */
    private final static String TOOL_GENE_MAPPING_SUBFOLDER = "resources/conf/gene_mappings/";
    /**
     * The name of the Ensembl versions file.
     */
    private static final String ENSEMBL_VERSIONS = "ensembl_versions";
    /**
     * The name of the GO domains file.
     */
    private static final String GO_DOMAINS = "go_domains";
    /**
     * The suffix to use for files containing gene mappings.
     */
    public final static String GENE_MAPPING_FILE_SUFFIX = "_gene_mappings";
    /**
     * The suffix to use for files containing GO mappings.
     */
    public final static String GO_MAPPING_FILE_SUFFIX = "_go_mappings";
    /**
     * The Ensembl versions for each species.
     */
    private HashMap<String, String> ensemblVersionsMap;
    /**
     * The horizontal padding to use when printing to the waiting dialog.
     */
    private final String PADDING = "    ";

    /**
     * Constructor.
     */
    public ProteinGeneDetailsProvider() {
    }

    /**
     * Initializes the factory. Note: the species factory must be initialized
     * first.
     *
     * @param jarFilePath the path to the jar file
     *
     * @throws java.io.IOException Exception thrown if an error occurs while
     * reading the species mapping
     */
    public void initialize(String jarFilePath) throws IOException {

        // load the previous ensembl version numbers
        File ensemblVersionsFile = getEnsemblVersionsFile();
        if (ensemblVersionsFile.exists()) {
            loadEnsemblSpeciesVersions(ensemblVersionsFile);
        } else {
            ensemblVersionsMap = new HashMap<>();
        }

        createDefaultGeneMappingFilesGeneric(
                jarFilePath,
                new File(jarFilePath, TOOL_GENE_MAPPING_SUBFOLDER + ENSEMBL_VERSIONS),
                new File(jarFilePath, TOOL_GENE_MAPPING_SUBFOLDER + GO_DOMAINS),
                true);
    }

    /**
     * Returns the gene maps for the given proteins. For every protein, the
     * species must be given as well as the gene name, in the format used in a
     * UniProt Importing gene mappings file.
     *
     * @param genePreferences the gene preferences
     * @param fastaSummary summary information on the Importing gene mappings
     * file containing the proteins
     * @param sequenceProvider the protein sequence provider
     * @param proteinDetailsProvider the protein details provider
     * @param waitingHandler waiting handler displaying progress for the
     * download and allowing canceling of the progress.
     *
     * @return the gene maps for the FASTA file loaded in the factory
     */
    public GeneMaps getGeneMaps(GeneParameters genePreferences, FastaSummary fastaSummary, SequenceProvider sequenceProvider, ProteinDetailsProvider proteinDetailsProvider, WaitingHandler waitingHandler) {

        Collection<String> accessions = sequenceProvider.getAccessions();

        SpeciesFactory speciesFactory = SpeciesFactory.getInstance();
        HashMap<String, GeneMapping> geneMappings = new HashMap<>(accessions.size());
        HashMap<String, GoMapping> goMappings = new HashMap<>(accessions.size());

        // download/update species mapping, put them in maps per species
        for (String uniprotTaxonomy : fastaSummary.speciesOccurrence.keySet()) {

            if (!uniprotTaxonomy.equals(SpeciesFactory.UNKNOWN)) {

                try {

                    Integer taxon = speciesFactory.getUniprotTaxonomy().getId(uniprotTaxonomy, true);

                    if (taxon != null) {

                        String speciesName = speciesFactory.getName(taxon);
                        String ensemblDatasetName = speciesFactory.getEnsemblDataset(taxon);

                        if (ensemblDatasetName != null) {

                            File geneMappingFile = getGeneMappingFile(ensemblDatasetName);
                            File goMappingFile = getGoMappingFile(ensemblDatasetName);

                            if (genePreferences.getAutoUpdate()) {
                                boolean success = true;
                                try {
                                    if (!geneMappingFile.exists() || !goMappingFile.exists() || newVersionExists(taxon)) { //@TODO: seems like the current version is never found
                                        success = downloadMappings(waitingHandler, taxon);
                                    }
                                    if (waitingHandler != null && waitingHandler.isRunCanceled()) {
                                        return null;
                                    }
                                } catch (Exception e) {
                                    e.printStackTrace();
                                    success = false;
                                }
                                if (!success) {
                                    waitingHandler.appendReport(PADDING + "Update of gene information for species " + speciesName + " failed.", true, true);
                                }
                            }

                            if (geneMappingFile.exists()) {
                                GeneMapping geneMapping = new GeneMapping();
                                try {
                                    geneMapping.importFromFile(geneMappingFile, waitingHandler);
                                    geneMappings.put(speciesName, geneMapping);
                                } catch (Exception e) {
                                    waitingHandler.appendReport(PADDING + "Import of gene mappings for " + speciesName + " failed.", true, true);
                                }
                            } else {
                                waitingHandler.appendReport(PADDING + "Gene mapping for " + speciesName + " not available.", true, true);
                            }

                            if (goMappingFile.exists()) {
                                GoMapping goMapping = new GoMapping();
                                try {
                                    goMapping.loadMappingsFromFile(goMappingFile, waitingHandler);
                                    goMappings.put(speciesName, goMapping);
                                } catch (Exception e) {
                                    waitingHandler.appendReport(PADDING + "Import of the GO mapping for " + speciesName + " failed.", true, true);
                                }
                            } else {
                                waitingHandler.appendReport(PADDING + "GO mapping for " + speciesName + " not available.", true, true);
                            }
                        } else {
                            waitingHandler.appendReport(PADDING + speciesName + " not available in Ensembl.", true, true);
                        }
                    }
                } catch (Exception e) {
                    waitingHandler.appendReport(PADDING + "No taxonomy found for " + uniprotTaxonomy + ".", true, true);
                }
            }
        }

        // get the mappings for the proteins in the sequence factory
        GeneMaps geneMaps = new GeneMaps();

        if (ensemblVersionsMap == null) {

            ensemblVersionsMap = new HashMap<>(1);

        }

        HashMap<String, String> ensemblVersionsUsed = new HashMap<>(ensemblVersionsMap);
        HashMap<String, String> geneNameToEnsemblIdMap = new HashMap<>(accessions.size());
        HashMap<String, String> geneNameToChromosomeMap = new HashMap<>(accessions.size());
        HashMap<String, HashSet<String>> proteinToGoMap = new HashMap<>(accessions.size());
        HashMap<String, HashSet<String>> goToProteinMap = new HashMap<>(accessions.size());
        HashMap<String, String> goNamesMap = new HashMap<>(accessions.size());

        for (String accession : accessions) {

            String uniprotTaxonomy = proteinDetailsProvider.getTaxonomy(accession);

            if (uniprotTaxonomy != null && !uniprotTaxonomy.equals("")) {

                try {
                    Integer taxon = speciesFactory.getUniprotTaxonomy().getId(uniprotTaxonomy, false);

                    if (taxon != null) {

                        String speciesName = speciesFactory.getName(taxon);

                        String geneName = proteinDetailsProvider.getGeneName(accession);

                        if (geneName != null) {

                            GeneMapping geneMapping = geneMappings.get(speciesName);

                            if (geneMapping != null) {

                                String chromosome = geneMapping.getChromosome(geneName);

                                if (chromosome != null) {

                                    geneNameToChromosomeMap.put(geneName, chromosome);

                                }

                                String ensemblId = geneMapping.getEnsemblAccession(geneName);

                                if (ensemblId != null) {

                                    geneNameToEnsemblIdMap.put(geneName, ensemblId);

                                }
                            }
                        }

                        GoMapping goMapping = goMappings.get(speciesName);

                        if (goMapping != null) {

                            HashSet<String> goTerms = proteinToGoMap.get(accession);

                            if (goTerms == null) {

                                goTerms = new HashSet<>(1);
                                proteinToGoMap.put(accession, goTerms);

                            }

                            HashSet<String> newTerms = goMapping.getGoAccessions(accession);

                            if (newTerms != null) {

                                goTerms.addAll(newTerms);

                                for (String goTerm : newTerms) {

                                    String goName = goMapping.getTermName(goTerm);

                                    if (goName != null) {

                                        goNamesMap.put(goTerm, goName);

                                    }

                                    HashSet<String> proteins = goToProteinMap.get(goTerm);

                                    if (proteins == null) {

                                        proteins = new HashSet<>(1);
                                        goToProteinMap.put(goTerm, proteins);

                                    }

                                    proteins.add(accession);

                                }
                            }
                        }
                    }
                } catch (Exception e) {

                    // Taxon not available, ignore
                    e.printStackTrace();

                }
            }
        }

        geneMaps.setEnsemblVersionsMap(ensemblVersionsUsed);
        geneMaps.setGeneNameToEnsemblIdMap(geneNameToEnsemblIdMap);
        geneMaps.setGeneNameToChromosomeMap(geneNameToChromosomeMap);
        geneMaps.setProteinToGoMap(proteinToGoMap);
        geneMaps.setGoAccessionToProteinMap(goToProteinMap);
        geneMaps.setGoNamesMap(goNamesMap);

        return geneMaps;
    }

    /**
     * Download the gene sequences mappings.
     *
     * @param destinationFile the destination file where to save the gene
     * sequences
     * @param ensemblType the Ensembl type, e.g., default or plants
     * @param ensemblSchemaName the Ensembl schema name, e.g., default or
     * plants_mart_18
     * @param ensemblDbName the Ensembl DB name of the selected species
     * @param waitingHandler waiting handler displaying progress and allowing
     * canceling the process
     *
     * @return true if downloading went OK
     *
     * @throws MalformedURLException if an MalformedURLException occurs
     * @throws IOException if an IOException occurs
     */
    public boolean downloadGeneSequences(File destinationFile, String ensemblType, String ensemblSchemaName, String ensemblDbName, WaitingHandler waitingHandler) throws MalformedURLException, IOException {

        // construct data
        String requestXml = "query=<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
                + "<!DOCTYPE Query>"
                + "<Query  virtualSchemaName = \"" + ensemblSchemaName + "\" formatter = \"FASTA\" header = \"0\" uniqueRows = \"1\" count = \"\" datasetConfigVersion = \"0.7\" >"
                + "<Dataset name = \"" + ensemblDbName + "\" interface = \"default\" >"
                + "<Attribute name = \"ensembl_gene_id\" />\n"
                + "<Attribute name = \"coding\" />"
                + "</Dataset>\n"
                + "</Query>"
                + "</Query>";

        String waitingText = "Downloading gene sequences. Please Wait...";

        return queryEnsembl(requestXml, waitingText, destinationFile, ensemblType, waitingHandler);
    }

    /**
     * Download the GO mappings.
     *
     * @param ensemblType the Ensembl type, e.g., default or plants
     * @param ensemblSchemaName the Ensembl schema name, e.g., default or
     * plants_mart_18
     * @param ensemblDbName the Ensembl db name of the selected species
     * @param swissProtMapping if true, use the uniprotswissprot_accession
     * parameter, if false use the uniprotsptrembl parameter
     * @param waitingHandler waiting handler displaying progress and allowing
     * canceling the process
     *
     * @return true if downloading went OK
     *
     * @throws MalformedURLException if an MalformedURLException occurs
     * @throws IOException if an IOException occurs
     */
    public boolean downloadGoMappings(String ensemblType, String ensemblSchemaName, String ensemblDbName, boolean swissProtMapping, WaitingHandler waitingHandler) throws MalformedURLException, IOException {

        String accessionMapping;

        if (swissProtMapping) {
            accessionMapping = "\"uniprotswissprot\"";
        } else {
            accessionMapping = "\"uniprotsptrembl\"";
        }

        // construct data
        String requestXml = "query=<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
                + "<!DOCTYPE Query>"
                + "<Query  virtualSchemaName = \"" + ensemblSchemaName + "\" formatter = \"TSV\" header = \"0\" uniqueRows = \"1\" count = \"\" datasetConfigVersion = \"0.7\" >"
                + "<Dataset name = \"" + ensemblDbName + "\" interface = \"default\" >"
                + "<Attribute name = " + accessionMapping + " />";

        requestXml += "<Attribute name = \"goslim_goa_accession\" />"
                + "<Attribute name = \"goslim_goa_description\" />";

        requestXml += "</Dataset>"
                + "</Query>";

        // @TODO: have to check if goslim_goa_accession and goslim_goa_description is available
        File tempFile = getGoMappingFile(ensemblDbName);

        String waitingText = "Downloading GO Mappings. Please Wait...";
        return queryEnsembl(requestXml, waitingText, tempFile, ensemblType, waitingHandler);
    }

    /**
     * Sends an XML query to Ensembl and writes the result in a text file.
     *
     * @param requestXml the XML request
     * @param destinationFile the file where to save the results
     * @param ensemblType the Ensembl type, e.g., default or plants
     *
     * @return true if downloading went OK
     *
     * @throws MalformedURLException if an MalformedURLException occurs
     * @throws IOException if an IOException occurs
     */
    public boolean queryEnsembl(String requestXml, File destinationFile, String ensemblType) throws MalformedURLException, IOException {
        return queryEnsembl(requestXml, destinationFile, ensemblType, null);
    }

    /**
     * Sends an XML query to Ensembl and writes the result in a text file.
     *
     * @param requestXml the XML request
     * @param destinationFile the file where to save the results
     * @param ensemblType the Ensembl type, e.g., default or plants
     * @param waitingHandler waiting handler displaying progress and allowing
     * canceling the process
     *
     * @return true if downloading went OK
     *
     * @throws MalformedURLException if an MalformedURLException occurs
     * @throws IOException if an IOException occurs
     */
    public boolean queryEnsembl(String requestXml, File destinationFile, String ensemblType, WaitingHandler waitingHandler) throws MalformedURLException, IOException {
        return queryEnsembl(requestXml, null, destinationFile, ensemblType, waitingHandler);
    }

    /**
     * Sends an XML query to Ensembl and writes the result in a text file.
     *
     * @param requestXml the XML request
     * @param destinationFile the file where to save the results
     * @param ensemblType the Ensembl type, e.g., default or plants
     * @param waitingHandler waiting handler displaying progress and allowing
     * canceling the process
     * @param waitingText the text to write in case a progress dialog is used
     *
     * @return true if downloading went OK
     *
     * @throws MalformedURLException if an MalformedURLException occurs
     * @throws IOException if an IOException occurs
     */
    public boolean queryEnsembl(String requestXml, String waitingText, File destinationFile, String ensemblType, WaitingHandler waitingHandler) throws MalformedURLException, IOException {

        if (waitingHandler != null && waitingHandler instanceof ProgressDialogX && waitingText == null) {
            waitingText = "Downloading from Ensembl. Please wait...";
        }
        boolean success = true;

        int lastThousand = 0;

        if (waitingHandler == null || !waitingHandler.isRunCanceled()) {

            // send data
            URL url = getEnsemblUrl(ensemblType);

            URLConnection conn = url.openConnection();
            conn.setDoOutput(true);
            OutputStreamWriter wr = new OutputStreamWriter(conn.getOutputStream());
            String lineBreak = System.getProperty("line.separator");

            try {
                wr.write(requestXml);
                wr.flush();

                if (waitingHandler == null || !waitingHandler.isRunCanceled()) {

                    if (waitingHandler != null) {
                        waitingHandler.setWaitingText(waitingText);
                    } else {
                        System.out.println(waitingText);
                    }

                    int counter = 0;

                    boolean fileCreated = destinationFile.createNewFile();

                    if (fileCreated || destinationFile.exists()) {

                        // get the response
                        BufferedReader br = new BufferedReader(new InputStreamReader(conn.getInputStream()));
                        try {
                            FileWriter w = new FileWriter(destinationFile);
                            try {
                                BufferedWriter bw = new BufferedWriter(w);

                                try {
                                    String rowLine = br.readLine();

                                    if (rowLine != null && rowLine.startsWith("Query ERROR")) {
                                        if (rowLine.lastIndexOf("Attribute goslim_goa_accession NOT FOUND") != -1) {
                                            success = false;
                                        } else if (rowLine.lastIndexOf("Attribute uniprotswissprot_accession NOT FOUND") != -1) {
                                            success = false;
                                        } else {
                                            throw new IllegalArgumentException("Query error: " + rowLine);
                                        }
                                    } else {
                                        while (rowLine != null && success) {
                                            if (waitingHandler != null) {
                                                if (waitingHandler.isRunCanceled()) {
                                                    break;
                                                }
                                                if (waitingHandler instanceof ProgressDialogX) {
                                                    waitingHandler.setWaitingText(waitingText + " (" + counter++ + " rows downloaded)");
                                                }
                                            } else {
                                                int thousand = ++counter / 10000;
                                                if (thousand > lastThousand) {
                                                    System.out.println(waitingText + " (" + counter + " rows downloaded)");
                                                    lastThousand = thousand;
                                                }
                                            }
                                            bw.write(rowLine + lineBreak);
                                            rowLine = br.readLine();
                                        }
                                    }
                                } finally {
                                    bw.close();
                                }
                            } finally {
                                w.close();
                            }
                        } finally {
                            br.close();
                        }
                    } else {
                        if (waitingHandler != null) {
                            waitingHandler.setRunCanceled();
                        }
                        throw new IllegalArgumentException("The mapping file could not be created.");
                    }
                }
            } finally {
                wr.close();
            }
        }

        return success;
    }

    /**
     * Download the gene mappings.
     *
     * @param ensemblType the Ensembl type, e.g., default or plants
     * @param ensemblSchemaName the Ensembl schema name, e.g., default or
     * plants_mart_18
     * @param ensemblDatasetName the Ensembl dataset name of the selected
     * species
     * @param ensemblVersion the Ensembl version
     * @param waitingHandler the waiting handler
     *
     * @throws MalformedURLException if an MalformedURLException occurs
     * @throws IOException if an IOException occurs
     * @throws IllegalArgumentException if an IllegalArgumentException occurs
     */
    public void downloadGeneMappings(String ensemblType, String ensemblSchemaName, String ensemblDatasetName, String ensemblVersion,
            WaitingHandler waitingHandler) throws MalformedURLException, IOException, IllegalArgumentException {

        // construct data
        String requestXml = "query=<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
                + "<!DOCTYPE Query>"
                + "<Query  virtualSchemaName = \"" + ensemblSchemaName + "\" formatter = \"TSV\" header = \"0\" uniqueRows = \"1\" count = \"\" datasetConfigVersion = \"0.7\" >"
                + "<Dataset name = \"" + ensemblDatasetName + "\" interface = \"default\" >"
                + "<Attribute name = \"ensembl_gene_id\" />"
                + "<Attribute name = \"external_gene_name\" />"
                + "<Attribute name = \"chromosome_name\" />"
                + "</Dataset>"
                + "</Query>";

        // @TODO: use the queryEnsembl method here as well?
        if (!waitingHandler.isRunCanceled()) {

            // send data
            URL url = getEnsemblUrl(ensemblType);

            URLConnection conn = url.openConnection();
            conn.setDoOutput(true);
            OutputStreamWriter wr = new OutputStreamWriter(conn.getOutputStream());
            String lineBreak = System.getProperty("line.separator");

            try {
                wr.write(requestXml);
                wr.flush();

                if (!waitingHandler.isRunCanceled()) {

                    waitingHandler.setWaitingText("Downloading Gene Mappings. Please Wait...");

                    int counter = 0;

                    File tempFile = getGeneMappingFile(ensemblDatasetName);

                    // get the response
                    BufferedReader br = new BufferedReader(new InputStreamReader(conn.getInputStream()));
                    try {
                        FileWriter w = new FileWriter(tempFile);
                        try {
                            BufferedWriter bw = new BufferedWriter(w);
                            try {
                                String rowLine = br.readLine();

                                if (rowLine != null && rowLine.startsWith("Query ERROR")) {
                                    throw new IllegalArgumentException("Query error on line: " + rowLine);
                                } else {
                                    while (rowLine != null && !waitingHandler.isRunCanceled()) {
                                        if (waitingHandler instanceof ProgressDialogX) {
                                            waitingHandler.setWaitingText("Downloading Gene Mappings. Please Wait... (" + counter++ + " rows downloaded)");
                                        }
                                        bw.write(rowLine + lineBreak);
                                        rowLine = br.readLine();
                                    }
                                }
                            } finally {
                                bw.close();
                            }
                        } finally {
                            w.close();
                        }
                    } finally {
                        br.close();
                    }

                    if (!waitingHandler.isRunCanceled()) {
                        updateEnsemblVersion(ensemblDatasetName, "Ensembl " + ensemblVersion);
                    }
                }
            } finally {
                wr.close();
            }
        }
    }

    /**
     * Returns the path to the folder containing the gene mapping files.
     *
     * @return the gene mapping folder
     */
    public static File getGeneMappingFolder() {
        return new File(GENE_MAPPING_FOLDER);
    }

    /**
     * Sets the folder where gene mappings are saved.
     *
     * @param geneMappingFolder the folder where gene mappings are saved
     */
    public static void setGeneMappingFolder(String geneMappingFolder) {
        GENE_MAPPING_FOLDER = geneMappingFolder;
    }

    /**
     * Copies the given gene mapping files to the gene mappings folder. If newer
     * versions of the mapping exists they will be overwritten according to
     * updateEqualVersion.
     *
     * @param jarFilePath the Ensembl versions file
     * @param sourceEnsemblVersionsFile the Ensembl versions file
     * @param sourceGoDomainsFile the GO domains file
     * @param updateEqualVersion if true, the version is updated with equal
     * version numbers, false, only update if the new version is newer
     */
    public void createDefaultGeneMappingFilesGeneric(String jarFilePath, File sourceEnsemblVersionsFile, File sourceGoDomainsFile, boolean updateEqualVersion) {

        if (!getGeneMappingFolder().exists()) {
            boolean folderCreated = getGeneMappingFolder().mkdirs();

            if (!folderCreated) {
                throw new IllegalArgumentException("Could not create the gene mapping folder.");
            }
        }

        File targetEnsemblVersionsFile = getEnsemblVersionsFile();
        File targetGoDomainsFile = getGoDomainsFile();

        HashMap<String, String> localEnsemblVersionsMap = new HashMap<String, String>();
        HashMap<String, Boolean> localUpdateSpeciesEnsembl = new HashMap<String, Boolean>();

        try {
            if (!targetEnsemblVersionsFile.exists()) {

                boolean fileCreated = targetEnsemblVersionsFile.createNewFile();

                if (!fileCreated) {
                    throw new IllegalArgumentException("Could not create the Ensembl versions file.");
                }

                IoUtil.copyFile(sourceEnsemblVersionsFile, targetEnsemblVersionsFile);

                localEnsemblVersionsMap = getEnsemblSpeciesVersions(targetEnsemblVersionsFile);
                for (Map.Entry<String, String> entry : localEnsemblVersionsMap.entrySet()) {
                    Integer speciesEnsemblVersionNew = getEnsemblVersionFromFile(sourceEnsemblVersionsFile, entry.getKey());
                    updateEnsemblVersion(entry.getKey(), "Ensembl " + speciesEnsemblVersionNew); // we actually just need to update the map
                    localUpdateSpeciesEnsembl.put(entry.getKey(), true);
                }
            } else {

                // file exists, just update every species ensembl version
                // read the "new" species Ensembl versions number
                localEnsemblVersionsMap = getEnsemblSpeciesVersions(targetEnsemblVersionsFile);

                for (Map.Entry<String, String> entry : localEnsemblVersionsMap.entrySet()) { // entry.getKey(): species ; entry.getValue() : version

                    Integer speciesEnsemblVersionNew = getEnsemblVersionFromFile(sourceEnsemblVersionsFile, entry.getKey());

                    if (speciesEnsemblVersionNew != null) {

                        Integer speciesEnsemblVersionOld = getEnsemblVersionFromFile(targetEnsemblVersionsFile, entry.getKey());

                        if (speciesEnsemblVersionOld == null
                                || speciesEnsemblVersionOld.equals(speciesEnsemblVersionNew) && updateEqualVersion
                                || speciesEnsemblVersionOld < speciesEnsemblVersionNew) {
                            localUpdateSpeciesEnsembl.put(entry.getKey(), true);
                            updateEnsemblVersion(entry.getKey(), "Ensembl " + speciesEnsemblVersionNew);
                        } else {
                            localUpdateSpeciesEnsembl.put(entry.getKey(), false);
                        }
                    } else {
                        localUpdateSpeciesEnsembl.put(entry.getKey(), false);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Could not create or update the Ensembl versions file.");
        }

        try {
            if (!targetGoDomainsFile.exists()) {
                boolean fileCreated = targetGoDomainsFile.createNewFile();
                if (!fileCreated) {
                    throw new IllegalArgumentException("Could not create the GO domains file.");
                }
            }
            IoUtil.copyFile(sourceGoDomainsFile, targetGoDomainsFile);
        } catch (IOException e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Could not create the GO domains file.");
        }

        for (Map.Entry<String, Boolean> entry : localUpdateSpeciesEnsembl.entrySet()) {
            if (entry.getValue()) {

                File sourceSpeciesGoMappingsFile = new File(jarFilePath, TOOL_GENE_MAPPING_SUBFOLDER + entry.getKey() + GO_MAPPING_FILE_SUFFIX);
                File sourceSpeciesGeneMappingFile = new File(jarFilePath, TOOL_GENE_MAPPING_SUBFOLDER + entry.getKey() + GENE_MAPPING_FILE_SUFFIX);
                File targetSpeciesGoMappingsFile = new File(getGeneMappingFolder(), sourceSpeciesGoMappingsFile.getName());
                File targetSpeciesGeneMappingFile = new File(getGeneMappingFolder(), sourceSpeciesGeneMappingFile.getName());

                try {

                    if (!targetSpeciesGoMappingsFile.exists()) {
                        boolean fileCreated = targetSpeciesGoMappingsFile.createNewFile();
                        if (!fileCreated) {
                            throw new IllegalArgumentException("Could not create the default species GO mapping file.");
                        }
                    }

                    IoUtil.copyFile(sourceSpeciesGoMappingsFile, targetSpeciesGoMappingsFile);

                } catch (IOException e) {
                    e.printStackTrace();
                    throw new IllegalArgumentException("Could not create the default species GO mapping file.");
                }

                try {
                    if (!targetSpeciesGeneMappingFile.exists()) {
                        boolean fileCreated = targetSpeciesGeneMappingFile.createNewFile();
                        if (!fileCreated) {
                            throw new IllegalArgumentException("Could not create the default species gene mapping file.");
                        }
                    }

                    IoUtil.copyFile(sourceSpeciesGeneMappingFile, targetSpeciesGeneMappingFile);

                } catch (IOException e) {
                    e.printStackTrace();
                    throw new IllegalArgumentException("Could not create the default species gene mapping file.");
                }
            }
        }
    }

    /**
     * Update the Ensembl version for the given species in the local map and in
     * the Ensembl versions file.
     *
     * @param ensemblDatasetName the dataset name of the species to update,
     * e.g., hsapiens_gene_ensembl
     * @param ensemblVersion the new Ensembl version
     *
     * @throws IOException if an IOException occurs
     */
    public void updateEnsemblVersion(String ensemblDatasetName, String ensemblVersion) throws IOException {

        if (ensemblVersionsMap == null) {
            ensemblVersionsMap = new HashMap<>();
        }
        ensemblVersionsMap.put(ensemblDatasetName, ensemblVersion);

        FileWriter w = new FileWriter(getEnsemblVersionsFile());
        try {
            BufferedWriter bw = new BufferedWriter(w);
            try {
                for (String key : ensemblVersionsMap.keySet()) {
                    bw.write(key + SEPARATOR + ensemblVersionsMap.get(key));
                    bw.newLine();
                }
            } finally {
                bw.close();
            }
        } finally {
            w.close();
        }
    }

    /**
     * Gets the Ensembl version of a given species from a file.
     *
     * @param ensemblVersionsFile the Ensembl versions file
     * @param species the species of interest
     *
     * @return the Ensembl version
     *
     * @throws IOException thrown whenever an error occurred while reading the
     * file
     */
    public Integer getEnsemblVersionFromFile(File ensemblVersionsFile, String species) throws IOException {
        Integer version = null;
        FileReader r = new FileReader(ensemblVersionsFile);
        try {
            BufferedReader br = new BufferedReader(r);
            try {
                String line;
                while ((line = br.readLine()) != null) {
                    String[] splittedLine = line.split(SEPARATOR);
                    String speciesAtLine = splittedLine[0];
                    if (speciesAtLine.equals(species)) {
                        String[] ensemblVersionSplit = splittedLine[1].split(" ");
                        version = new Integer(ensemblVersionSplit[1]);
                    }
                }
            } finally {
                br.close();
            }
        } finally {
            r.close();
        }
        return version;
    }

    /**
     * Gets the information contained into the Ensembl species file.
     * 
     * @param ensemblVersionsFile the Ensembl species file to read
     * @return The Ensembl versions for each species
     * @throws FileNotFoundException if an FileNotFoundException occurs
     * @throws IOException if an IOException occurs
     */
    public HashMap<String, String> getEnsemblSpeciesVersions(File ensemblVersionsFile) throws FileNotFoundException, IOException {

        HashMap<String, String> localEnsemblVersionsMap = new HashMap<String, String>();

        // load the existing ensembl version numbers
        FileReader r = new FileReader(ensemblVersionsFile);

        try {
            BufferedReader br = new BufferedReader(r);

            try {
                String line = br.readLine();

                while (line != null) {
                    String[] elements = line.split("\\t");
                    localEnsemblVersionsMap.put(elements[0], elements[1]);
                    line = br.readLine();
                }
            } finally {
                br.close();
            }
        } finally {
            r.close();
        }

        return localEnsemblVersionsMap;
    }

    /**
     * Loads the given Ensembl species file.
     *
     * @param ensemblVersionsFile the Ensembl species file to load
     *
     * @throws IOException if an IOException occurs
     */
    public void loadEnsemblSpeciesVersions(File ensemblVersionsFile) throws IOException {

        // load the existing ensembl version numbers
        try (BufferedReader br = new BufferedReader(new FileReader(ensemblVersionsFile))) {

            ensemblVersionsMap = new HashMap<>();
            String line = br.readLine();

            while (line != null) {
                String[] elements = line.split("\\t");
                ensemblVersionsMap.put(elements[0], elements[1]);
                line = br.readLine();
            }
        }
    }

    /**
     * Returns the Ensembl URL for the given Ensembl (sub-)version.
     *
     * @param ensemblType the Ensembl type, e.g., fungi or plants
     *
     * @return the Ensembl URL
     *
     * @throws MalformedURLException exception thrown if the URL is malformed
     */
    private URL getEnsemblUrl(String ensemblType) throws MalformedURLException {

        if (ensemblType.equalsIgnoreCase("fungi")) {

            return new URL("https://fungi.ensembl.org/biomart/martservice/result");

        } else if (ensemblType.equalsIgnoreCase("plants")) {

            return new URL("https://plants.ensembl.org/biomart/martservice/result");

        } else if (ensemblType.equalsIgnoreCase("protists")) {

            return new URL("https://protists.ensembl.org/biomart/martservice/result");

        } else if (ensemblType.equalsIgnoreCase("metazoa")) {

            return new URL("https://metazoa.ensembl.org/biomart/martservice/result");

        } else {

            return new URL("https://www.ensembl.org/biomart/martservice/result");

        }
    }

    /**
     * Try to download the gene and GO mappings for the currently selected
     * species.
     *
     * @param waitingHandler the waiting handler
     * @param taxon the NCBI taxon of the species
     *
     * @return true if the download was successful
     *
     * @throws java.io.IOException exception thrown whenever an error occurred
     * while reading the mapping files
     */
    public boolean downloadMappings(WaitingHandler waitingHandler, Integer taxon) throws IOException {

        SpeciesFactory speciesFactory = SpeciesFactory.getInstance();

        String latinName = speciesFactory.getLatinName(taxon);
        if (latinName == null) {
            latinName = taxon.toString();
        }

        if (waitingHandler.isReport()) {
            waitingHandler.appendReport(PADDING + "Downloading GO and gene mappings for species " + latinName + ".", true, true);
        }

        EnsemblGenomeDivision ensemblGenomeDivision = speciesFactory.getEnsemblGenomesSpecies().getDivision(taxon);

        String ensemblType = "ensembl";
        if (ensemblGenomeDivision != null) {
            ensemblType = ensemblGenomeDivision.ensemblType;
        }

        String schemaName = EnsemblVersion.getEnsemblSchemaName(ensemblGenomeDivision);
        if (schemaName == null) {
            return false;
        }

        String ensemblDatasetName = speciesFactory.getEnsemblDataset(taxon);
        if (ensemblDatasetName == null) {
            return false;
        }

        if (!waitingHandler.isRunCanceled()) {

            boolean goMappingsDownloaded = downloadGoMappings(ensemblType, schemaName, ensemblDatasetName, true, waitingHandler);

            // swiss prot mapping not found, try trembl
            if (!goMappingsDownloaded) {
                goMappingsDownloaded = downloadGoMappings(ensemblType, schemaName, ensemblDatasetName, false, waitingHandler);
            }

            if (!goMappingsDownloaded) {
                waitingHandler.appendReport(PADDING + "Gene ontology mappings not available. Downloading gene mappings only.", true, true);
            } else if (waitingHandler.isReport()) {
                waitingHandler.appendReport(PADDING + "GO mappings downloaded.", true, true);
            }
        }

        if (!waitingHandler.isRunCanceled()) {
            downloadGeneMappings(ensemblType, schemaName, ensemblDatasetName,
                    EnsemblVersion.getCurrentEnsemblVersion(ensemblGenomeDivision).toString(), waitingHandler);

            if (!waitingHandler.isRunCanceled()) {
                if (waitingHandler.isReport()) {
                    waitingHandler.appendReport(PADDING + "Gene mappings downloaded.", true, true);
                }
            }
        }

        boolean canceled = waitingHandler.isRunCanceled();
        return !canceled;
    }

    /**
     * Returns the gene mapping file.
     *
     * @param ensemblDatasetName the Ensembl dataset name
     *
     * @return the gene mapping file
     */
    public static File getGeneMappingFile(String ensemblDatasetName) {
        return new File(getGeneMappingFolder(), ensemblDatasetName + GENE_MAPPING_FILE_SUFFIX);
    }

    /**
     * Returns the GO mapping file.
     *
     * @param ensemblDatasetName the Ensembl dataset name
     *
     * @return the GO mapping file
     */
    public static File getGoMappingFile(String ensemblDatasetName) {
        return new File(getGeneMappingFolder(), ensemblDatasetName + GO_MAPPING_FILE_SUFFIX);
    }

    /**
     * Returns the Ensembl version file.
     *
     * @return the Ensembl version file
     */
    public static File getEnsemblVersionsFile() {
        return new File(getGeneMappingFolder(), ENSEMBL_VERSIONS);
    }

    /**
     * Returns the GO domains file.
     *
     * @return the GO domains file
     */
    public static File getGoDomainsFile() {
        return new File(getGeneMappingFolder(), GO_DOMAINS);
    }

    /**
     * Returns the Ensembl version for a given species.
     *
     * @param taxon the NCBI taxon of the species
     *
     * @return the Ensembl version for a given species.
     */
    public String getEnsemblVersion(Integer taxon) {
        SpeciesFactory speciesFactory = SpeciesFactory.getInstance();
        String ensemblDatasetName = speciesFactory.getEnsemblDataset(taxon);
        if (ensemblVersionsMap == null) {
            return null;
        }
        return ensemblVersionsMap.get(ensemblDatasetName);
    }

    /**
     * Returns true if a newer version of the species mapping exists in Ensembl.
     *
     * @param taxon the NCBI taxon of the species
     *
     * @return rue if a newer version of the species mapping exists in Ensemble
     */
    public boolean newVersionExists(Integer taxon) {

        EnsemblGenomeDivision ensemblGenomeDivision = SpeciesFactory.getInstance().getEnsemblGenomesSpecies().getDivision(taxon);

        Integer latestEnsemblVersion = EnsemblVersion.getCurrentEnsemblVersion(ensemblGenomeDivision);
        String currentEnsemblVersionAsString = getEnsemblVersion(taxon);

        if (currentEnsemblVersionAsString != null) {

            currentEnsemblVersionAsString = currentEnsemblVersionAsString.substring(currentEnsemblVersionAsString.indexOf(" ") + 1);
            Integer currentEnsemblVersion;

            try {
                currentEnsemblVersion = new Integer(currentEnsemblVersionAsString);
            } catch (NumberFormatException e) {
                e.printStackTrace();
                currentEnsemblVersion = latestEnsemblVersion;
            }

            return currentEnsemblVersion < latestEnsemblVersion;
        }

        return true;
    }
}
