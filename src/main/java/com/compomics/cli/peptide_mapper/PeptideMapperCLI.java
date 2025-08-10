package com.compomics.cli.peptide_mapper;

import com.compomics.util.experiment.biology.aminoacids.sequence.AminoAcidSequence;
import com.compomics.util.experiment.identification.amino_acid_tags.MassGap;
import com.compomics.util.experiment.identification.amino_acid_tags.Tag;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.parameters.identification.search.ModificationParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.experiment.identification.protein_inference.fm_index.FMIndex;
import com.compomics.util.gui.waiting.waitinghandlers.WaitingHandlerCLIImpl;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import java.io.File;
import java.io.PrintWriter;
import com.compomics.util.experiment.identification.protein_inference.FastaMapper;
import com.compomics.util.experiment.io.identification.IdfileReader;
import com.compomics.util.experiment.io.identification.IdfileReaderFactory;
import com.compomics.util.experiment.io.mass_spectrometry.MsFileHandler;
import com.compomics.util.io.IoUtil;
import com.compomics.util.io.compression.ZipUtils;
import com.compomics.util.io.flat.SimpleFileReader;
import com.compomics.util.io.flat.SimpleFileWriter;
import com.compomics.util.threading.SimpleSemaphore;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.*;
import java.util.stream.Collectors;

/**
 * Command line peptide mapping.
 *
 * @author Dominik Kopczynski
 */
public class PeptideMapperCLI {

    public static int TIMEOUT_DAYS = 1;

    public static void printHelp() {
        System.out.println("PeptideMapperCLI: a tool to map peptides or sequence tags against a given proteome.");
        System.out.println("usage: PeptideMapperCLI -[p|t|c] input-fasta input-peptide/tag-csv output-csv [additonal options]");
        System.out.println();
        System.out.println("Options are:");
        System.out.println("\t-p\tpeptide mapping");
        System.out.println("\t-t\tsequence tag mapping");
        System.out.println("\t-x\textract peptide spectrum matches from SearchGUI output file");
        System.out.println("\t-h\tprint this info");
        System.out.println();
        System.out.println("Additional options are:");
        System.out.println("\t-u [utilities-parameter-file]\tpeptide mapping");
        System.out.println("\t-f\tadd flanking amino acids to peptide in output");
        System.out.println("\t-c\tspecify the number of cores used");

        System.out.println();
        System.out.println("Default parameters:");
        System.out.println("\tindexing method:\t\tfm-index");
        System.out.println("\tframentation tolerance [Da]:\t0.02");
    }

    /**
     * Main class.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {

        if (args.length < 1 || args[0].equals("-h") || args[0].equals("--help")) {
            printHelp();
            System.exit(-1);
        }

        switch (args[0]) {
            case "-x":
                convertSearchGUIFile(args);
                break;

            case "-p":
            case "-t":
                handleParameters(args);
                break;

            default:
                printHelp();
                System.exit(-1);

        }
    }

    public static String tagToString(Tag tag) {

        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < tag.getContent().size(); ++i) {
            if (i > 0) {
                sb.append(",");
            }
            if (tag.getContent().get(i) instanceof MassGap) {
                sb.append(tag.getContent().get(i).getMass());
            } else if (tag.getContent().get(i) instanceof AminoAcidSequence) {
                sb.append(tag.getContent().get(i).asSequence());
            }
        }

        return sb.toString();
    }

    public static void convertSearchGUIFile(String[] args) {

        if (args.length < 2) {
            printHelp();
            System.exit(-1);
        }

        IdentificationParameters identificationParameters = new IdentificationParameters();
        boolean parametersLoaded = false;

        ArrayList<File> idFiles = new ArrayList<>();

        try {
            File SGfile = new File(args[1]);
            MsFileHandler msFileHandler = new MsFileHandler();
            File tmpDir = Paths.get(SGfile.getAbsoluteFile().getParent(), "tmpDir").toFile();
            tmpDir.mkdir();
            ZipUtils.unzip(
                    SGfile,
                    tmpDir,
                    null
            );
            IdfileReaderFactory readerFactory = IdfileReaderFactory.getInstance();

            ArrayList<File> dirs = new ArrayList<>();
            dirs.add(tmpDir);

            for (int d = 0; d < dirs.size(); ++d) {
                File dir = dirs.get(d);
                for (File zippedFile : dir.listFiles()) {
                    if (zippedFile.isDirectory()) {
                        dirs.add(zippedFile);
                    }
                    String lowerCaseName = zippedFile.getName().toLowerCase();

                    if (lowerCaseName.endsWith(".dat")
                            || lowerCaseName.endsWith(".omx")
                            || lowerCaseName.endsWith(".res")
                            || lowerCaseName.endsWith(".xml")
                            || lowerCaseName.endsWith(".mzid")
                            || lowerCaseName.endsWith(".csv")
                            || lowerCaseName.endsWith(".tsv")
                            || lowerCaseName.endsWith(".tags")
                            || lowerCaseName.endsWith(".pnovo.txt")
                            || lowerCaseName.endsWith(".tide-search.target.txt")
                            || lowerCaseName.endsWith(".psm.gz")
                            || lowerCaseName.endsWith(".omx.gz")
                            || lowerCaseName.endsWith(".res.gz")
                            || lowerCaseName.endsWith(".xml.gz")
                            || lowerCaseName.endsWith(".mzid.gz")
                            || lowerCaseName.endsWith(".csv.gz")
                            || lowerCaseName.endsWith(".tsv.gz")
                            || lowerCaseName.endsWith(".tags.gz")
                            || lowerCaseName.endsWith(".pnovo.txt.gz")
                            || lowerCaseName.endsWith(".tide-search.target.txt.gz")
                            || lowerCaseName.endsWith(".psm.gz")) {

                        if (!lowerCaseName.endsWith("mods.xml")
                                && !lowerCaseName.endsWith("usermods.xml")
                                && !lowerCaseName.endsWith("settings.xml")) {

                            idFiles.add(zippedFile);

                        }
                    } else if (lowerCaseName.endsWith(".par")) {

                        try {
                            identificationParameters = IdentificationParameters.getIdentificationParameters(zippedFile);
                            parametersLoaded = true;
                        } catch (Exception e) {
                            System.err.println("Error: cound not open or parse parameter file");
                            System.exit(-1);
                        }

                    } else if (lowerCaseName.endsWith(".mgf")) {
                        msFileHandler.register(zippedFile, zippedFile.getParentFile(), null);
                    }
                }
            }

            if (!parametersLoaded) {
                identificationParameters.getSearchParameters().setModificationParameters(new ModificationParameters());
                identificationParameters.getSearchParameters().setFragmentIonAccuracy(0.02);
                identificationParameters.getSearchParameters().setFragmentAccuracyType(SearchParameters.MassAccuracyType.DA);

                identificationParameters.getSequenceMatchingParameters().setSequenceMatchingType(SequenceMatchingParameters.MatchingType.indistiguishableAminoAcids);
                identificationParameters.getSequenceMatchingParameters().setLimitX(0.25);

            }

            ArrayList<SpectrumMatch> spectrumMatches = new ArrayList<>();

            for (File idFile : idFiles) {

                IdfileReader fileReader = null;
                try {

                    fileReader = readerFactory.getFileReader(idFile);

                } catch (OutOfMemoryError error) {
                    System.out.println("Ran out of memory when parsing \'" + IoUtil.getFileName(idFile) + "\'.");
                    System.exit(-1);
                }

                spectrumMatches.addAll(fileReader.getAllSpectrumMatches(
                        msFileHandler,
                        null,
                        identificationParameters.getSearchParameters(),
                        identificationParameters.getSequenceMatchingParameters(),
                        true
                ));

            }

            boolean hasPeptides = false;
            boolean hasTags = false;
            for (SpectrumMatch spectrumMatch : spectrumMatches) {
                hasPeptides |= spectrumMatch.hasPeptideAssumption();
                hasTags |= spectrumMatch.hasTagAssumption();
            }

            PrintWriter writerPeptides = null;
            PrintWriter writerTags = null;

            try {
                if (hasPeptides) {
                    writerPeptides = new PrintWriter(Paths.get((new File(args[1])).getParent(), "extracted-peptides.csv").toFile(), "UTF-8");
                }
                if (hasTags) {
                    writerTags = new PrintWriter(Paths.get((new File(args[1])).getParent(), "extracted-tags.csv").toFile(), "UTF-8");
                }
            } catch (Exception e) {
                handleError("Opening error", "Error: could not open output files properly.", e);
            }

            for (SpectrumMatch spectrumMatch : spectrumMatches) {
                if (spectrumMatch.hasPeptideAssumption()) {
                    ArrayList<String> peptides = spectrumMatch.getAllPeptideAssumptions()
                            .map(assumption -> assumption.getPeptide().getSequence())
                            .collect(Collectors.toCollection(ArrayList::new));

                    for (String peptide : peptides) {
                        try {
                            writerPeptides.println(peptide);
                        } catch (Exception e) {
                            handleError("Writing error", "Error: could not write into file.", e);
                        }
                    }
                } else if (spectrumMatch.hasTagAssumption()) {
                    ArrayList<String> sequenceTags = spectrumMatch.getAllTagAssumptions()
                            .map(assumption -> tagToString(assumption.getTag()))
                            .collect(Collectors.toCollection(ArrayList::new));

                    for (String sequenceTag : sequenceTags) {
                        try {
                            writerTags.println(sequenceTag);
                        } catch (Exception e) {
                            handleError("Writing error", "Error: could not write into file.", e);
                        }
                    }
                }
            }

            // close everything
            try {
                if (writerPeptides != null) {
                    writerPeptides.close();
                }
                if (writerTags != null) {
                    writerTags.close();
                }
            } catch (Exception e) {
                handleError("Closing error", "Error: could not close files properly", e);
            }
            IoUtil.deleteDir(tmpDir);

        } catch (Exception e) {
            handleError("SearchGUI file error", "Could not open SearchGUI zip file properly.", e);
        }
    }

    public static void handleParameters(String[] args) {

        if (args.length < 4) {
            printHelp();
            System.exit(-1);
        }

        WaitingHandlerCLIImpl waitingHandlerCLIImpl = new WaitingHandlerCLIImpl();
        File fastaFile = new File(args[1]);
        String inputFileName = args[2];
        String outputFileName = args[3];
        boolean flanking = false;
        boolean peptideMapping = args[0].equals("-p");
        int nCores = Runtime.getRuntime().availableProcessors();

        if (!args[0].equals("-p") && !args[0].equals("-t")) {
            System.out.println("Invalid first parameter: " + args[0]);
            System.exit(-1);
        }

        // read in the parameters
        IdentificationParameters identificationParameters = new IdentificationParameters(new SearchParameters());
        boolean customParameters = false;
        if (args.length >= 5) {
            int argPos = 4;
            while (argPos < args.length) {
                switch (args[argPos]) {
                    case "-f":  // flanking
                        flanking = true;
                        ++argPos;
                        break;

                    case "-c": // number of cores
                        try {
                            nCores = Integer.parseInt(args[argPos + 1]);
                            if (nCores < 1 || 100000 < nCores) {
                                throw new Exception();
                            }
                        } catch (Exception e) {
                            System.out.println("Parameter -c has no valid number.");
                            System.exit(-1);
                        }
                        argPos += 2;
                        break;

                    case "-u":  // use utilities parameter file

                        try {
                            File parameterFile = new File(args[argPos + 1]);
                            identificationParameters = IdentificationParameters.getIdentificationParameters(parameterFile);
                        } catch (Exception e) {
                            System.err.println("Error: cound not open or parse parameter file");
                            System.exit(-1);
                        }
                        customParameters = true;
                        argPos += 2;
                        break;

                    default:
                        ++argPos;
                        break;
                }
            }

        }
        if (!customParameters) {
            identificationParameters.getSearchParameters().setModificationParameters(new ModificationParameters());
            identificationParameters.getSearchParameters().setFragmentIonAccuracy(0.02);
            identificationParameters.getSearchParameters().setFragmentAccuracyType(SearchParameters.MassAccuracyType.DA);
            identificationParameters.getSequenceMatchingParameters().setSequenceMatchingType(SequenceMatchingParameters.MatchingType.indistiguishableAminoAcids);
            identificationParameters.getSequenceMatchingParameters().setLimitX(0.25);
        }
        identificationParameters.getSearchParameters().setFlanking(flanking);

        runMapping(fastaFile,
                waitingHandlerCLIImpl,
                identificationParameters,
                inputFileName,
                outputFileName,
                nCores,
                peptideMapping);
    }

    public static void handleError(String outputFileName, String errorMessage, Throwable e) {
        PrintWriter writer = null;
        System.out.println(errorMessage);

        String path = (new File(outputFileName)).getParent();
        String timeStamp = new SimpleDateFormat("yyyy-MM-dd_HH-mm").format(new Date());

        String logFileName = Paths.get(path, "PeptideMapper_Error-log_" + timeStamp + ".txt").toString();

        try {
            if (e instanceof OutOfMemoryError) {
                System.gc();
            }
            writer = new PrintWriter(logFileName, "UTF-8");
            e.printStackTrace(writer);
            writer.close();
            System.out.println("A proper error log was stored in '" + logFileName + "'");

        } catch (Exception e2) {
            System.out.println("Could not store the error log into a file, it will be printed directly:");
            System.out.println();
            System.out.println(e.getMessage());
        }
        System.exit(-1);
    }

    public static void runMapping(File fastaFile,
            WaitingHandlerCLIImpl waitingHandlerCLIImpl,
            IdentificationParameters identificationParameters,
            String inputFilePath,
            String outputFilePath,
            int nCores,
            boolean peptideMapping) {

        // setting up the mapper
        FastaMapper peptideMapper = null;
        System.out.println("Start indexing fasta file");
        long startTimeIndex = System.nanoTime();
        try {
            peptideMapper = new FMIndex(fastaFile, null, waitingHandlerCLIImpl, true, identificationParameters);
        } catch (Exception e) {
            handleError(outputFilePath, "Error: cound not index the fasta file", e);
        } catch (OutOfMemoryError e) {
            handleError(outputFilePath, "Error: not enough memory available. Please try to run Java with more memory, e.g.: java -Xmx16G ...\n\n", e);
        }

        double diffTimeIndex = System.nanoTime() - startTimeIndex;
        System.out.println();
        System.out.println("Indexing took " + (diffTimeIndex / 1e9) + " seconds and consumes " + (((float) ((FMIndex) peptideMapper).getAllocatedBytes()) / 1e6) + " MB");
        System.out.println();
        System.out.println("Start mapping using " + nCores + " threads");

        // Set up progress
        long lineCount = 0;
        try {
            Path path = Paths.get(inputFilePath);
            lineCount = Files.lines(path).count();
        } catch (Exception e) {
            handleError(outputFilePath, "Error: could not open files properly.", e);
        }

        waitingHandlerCLIImpl.setSecondaryProgressCounterIndeterminate(false);
        waitingHandlerCLIImpl.setMaxSecondaryProgressCounter((int) lineCount);
        waitingHandlerCLIImpl.setSecondaryProgressCounter(0);

        // Set up reader and writer
        // starting the mapping
        try (SimpleFileReader simpleFileReader = SimpleFileReader.getFileReader(new File(inputFilePath))) {

            try (SimpleFileWriter writer = new SimpleFileWriter(new File(outputFilePath), false)) {

                long startTimeMapping = System.nanoTime();

                ArrayList<MappingWorker> workers = new ArrayList<>();
                ExecutorService importPool = Executors.newFixedThreadPool(nCores);
                SimpleSemaphore readMutex = new SimpleSemaphore(1);

                for (int i = 0; i < nCores; ++i) {
                    MappingWorker mw = new MappingWorker(
                            waitingHandlerCLIImpl,
                            peptideMapper,
                            identificationParameters,
                            simpleFileReader,
                            readMutex,
                            writer,
                            peptideMapping
                    );
                    importPool.submit(mw);
                    workers.add(mw);
                }
                importPool.shutdown();
                if (!importPool.awaitTermination(TIMEOUT_DAYS, TimeUnit.DAYS)) {
                    System.out.println("Analysis timed out (time out: " + TIMEOUT_DAYS + " days)");
                }
                for (MappingWorker mw : workers) {
                    if (mw.exception != null) {
                        throw mw.exception;
                    }
                }

                long diffTimeMapping = System.nanoTime() - startTimeMapping;
                System.out.println();
                System.out.println("Mapping " + lineCount + " peptides took " + (diffTimeMapping / 1e9) + " seconds");

            }
        } catch (Exception e) {
            handleError(outputFilePath, "Error: mapping went wrong", e);
        }
    }
}
