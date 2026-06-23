package com.compomics.util.experiment.io.identification.idfilereaders;

import com.compomics.util.Util;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.identification.spectrum_assumptions.PeptideAssumption;
import com.compomics.util.experiment.io.identification.IdfileReader;
import com.compomics.util.experiment.mass_spectrometry.SpectrumProvider;
import com.compomics.util.io.IoUtil;
import com.compomics.util.io.flat.SimpleFileReader;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.waiting.WaitingHandler;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.xml.bind.JAXBException;

/**
 * Shared parser for InstaNovo v1.2.2 normalized CSV predictions.
 *
 * @author CompOmics
 */
abstract class InstaNovoCsvIdfileReader implements IdfileReader {

    /**
     * The supported InstaNovo version.
     */
    private static final String SOFTWARE_VERSION = "1.2.2";
    /**
     * Pattern matching common scan or index tokens in spectrum titles.
     */
    private static final Pattern TITLE_NUMBER_PATTERN = Pattern.compile("(?i)(?:scan|index|scan_number)\\s*[=: ]\\s*(\\d+)");
    /**
     * The CSV file.
     */
    private final File csvFile;
    /**
     * The advocate used for peptide assumptions.
     */
    private final Advocate advocate;
    /**
     * The extension this reader is registered for.
     */
    private final String extension;

    /**
     * Constructor.
     *
     * @param csvFile the CSV file
     * @param advocate the advocate
     * @param extension the registered extension
     */
    protected InstaNovoCsvIdfileReader(File csvFile, Advocate advocate, String extension) {
        this.csvFile = csvFile;
        this.advocate = advocate;
        this.extension = extension;
    }

    @Override
    public String getExtension() {
        return extension;
    }

    @Override
    public ArrayList<SpectrumMatch> getAllSpectrumMatches(
            SpectrumProvider spectrumProvider,
            WaitingHandler waitingHandler,
            SearchParameters searchParameters
    )
            throws IOException, IllegalArgumentException, SQLException, ClassNotFoundException, InterruptedException, JAXBException {

        return getAllSpectrumMatches(
                spectrumProvider,
                waitingHandler,
                searchParameters,
                null,
                true
        );
    }

    @Override
    public ArrayList<SpectrumMatch> getAllSpectrumMatches(
            SpectrumProvider spectrumProvider,
            WaitingHandler waitingHandler,
            SearchParameters searchParameters,
            SequenceMatchingParameters sequenceMatchingPreferences,
            boolean expandAaCombinations
    )
            throws IOException, IllegalArgumentException, SQLException, ClassNotFoundException, InterruptedException, JAXBException {

        if (spectrumProvider == null) {
            throw new IllegalArgumentException("A spectrum provider is required to import InstaNovo results.");
        }

        ArrayList<SpectrumMatch> result = new ArrayList<>();
        HashMap<String, SpectrumMatch> matches = new HashMap<>();
        HashMap<String, SpectrumTitleLookup> spectrumTitleLookups = new HashMap<>();

        try (SimpleFileReader reader = SimpleFileReader.getFileReader(csvFile)) {

            String line = reader.readLine();

            if (line == null) {
                throw new IllegalArgumentException("The InstaNovo csv file is empty.");
            }

            ArrayList<String> headers = parseCsvLine(line);
            HashMap<String, Integer> columnIndexes = getColumnIndexes(headers);

            int experimentIndex = getOptionalColumn(columnIndexes, "experiment_name");
            int spectrumIdIndex = getOptionalColumn(columnIndexes, "spectrum_id", "spectrum");
            int scanNumberIndex = getOptionalColumn(columnIndexes, "scan_number", "scan");
            int chargeIndex = getRequiredColumn(columnIndexes, "precursor_charge", "charge", "z");
            int predictionIndex = getRequiredColumn(columnIndexes, "predictions", "prediction", "sequence");
            int scoreIndex = getRequiredColumn(columnIndexes, "log_probs", "prediction_log_probability", "predictions_log_probability");

            if (experimentIndex < 0 && spectrumIdIndex < 0 && scanNumberIndex < 0) {
                throw new IllegalArgumentException("Mandatory spectrum identification columns are missing in the InstaNovo csv file.");
            }

            int lineNumber = 1;
            while ((line = reader.readLine()) != null) {

                lineNumber++;

                if (line.trim().isEmpty()) {
                    continue;
                }

                ArrayList<String> values = parseCsvLine(line);

                String prediction = getValue(values, predictionIndex).trim();

                if (prediction.isEmpty()) {
                    continue;
                }

                String experimentName = experimentIndex >= 0 ? getValue(values, experimentIndex).trim() : "";
                String spectrumId = spectrumIdIndex >= 0 ? getValue(values, spectrumIdIndex).trim() : "";
                String scanNumber = scanNumberIndex >= 0 ? getValue(values, scanNumberIndex).trim() : "";
                Integer charge = getCharge(getValue(values, chargeIndex), lineNumber, waitingHandler);

                if (charge == null) {
                    continue;
                }

                String spectrumFileName = getSpectrumFileName(spectrumProvider, experimentName, spectrumId);
                SpectrumTitleLookup spectrumTitleLookup = spectrumTitleLookups.get(spectrumFileName);

                if (spectrumTitleLookup == null) {
                    spectrumTitleLookup = new SpectrumTitleLookup(spectrumProvider, spectrumFileName);
                    spectrumTitleLookups.put(spectrumFileName, spectrumTitleLookup);
                }

                String spectrumTitle = getSpectrumTitle(spectrumTitleLookup, spectrumFileName, spectrumId, scanNumber);

                double logProbability = Util.readDoubleAsString(getValue(values, scoreIndex));
                double score = logProbability;

                ParsedPeptide parsedPeptide = parsePeptide(prediction, lineNumber);
                Peptide peptide = new Peptide(parsedPeptide.sequence, parsedPeptide.modificationMatches);
                peptide.estimateTheoreticMass(
                        searchParameters.getModificationParameters(),
                        null,
                        SequenceMatchingParameters.DEFAULT_STRING_MATCHING
                );
                PeptideAssumption peptideAssumption = new PeptideAssumption(
                        peptide,
                        1,
                        advocate.getIndex(),
                        charge,
                        logProbability,
                        score,
                        IoUtil.getFileName(csvFile)
                );

                String matchKey = spectrumFileName + "\n" + spectrumTitle;
                SpectrumMatch spectrumMatch = matches.get(matchKey);

                if (spectrumMatch == null) {
                    spectrumMatch = new SpectrumMatch(spectrumFileName, spectrumTitle);
                    matches.put(matchKey, spectrumMatch);
                    result.add(spectrumMatch);
                }

                spectrumMatch.addPeptideAssumption(advocate.getIndex(), peptideAssumption);
            }
        }

        return result;
    }

    @Override
    public void close() throws IOException {
        // Nothing to close.
    }

    @Override
    public HashMap<String, ArrayList<String>> getSoftwareVersions() {

        HashMap<String, ArrayList<String>> result = new HashMap<>();
        ArrayList<String> versions = new ArrayList<>();
        versions.add(SOFTWARE_VERSION);
        result.put(advocate.getName(), versions);

        if (advocate == Advocate.instanovoRefined) {

            ArrayList<String> instaNovoVersions = new ArrayList<>();
            instaNovoVersions.add(SOFTWARE_VERSION);
            result.put(Advocate.instanovo.getName(), instaNovoVersions);

            ArrayList<String> instaNovoPlusVersions = new ArrayList<>();
            instaNovoPlusVersions.add(SOFTWARE_VERSION);
            result.put(Advocate.instanovoPlus.getName(), instaNovoPlusVersions);
        }

        return result;
    }

    @Override
    public boolean hasDeNovoTags() {
        return false;
    }

    /**
     * Returns the spectrum file name without extension.
     *
     * @param spectrumProvider the spectrum provider
     * @param experimentName the experiment name
     * @param spectrumId the spectrum id
     *
     * @return the spectrum file name without extension
     */
    private String getSpectrumFileName(SpectrumProvider spectrumProvider, String experimentName, String spectrumId) {

        String fileName = experimentName;

        if (fileName == null || fileName.isEmpty()) {
            int separatorIndex = spectrumId.indexOf(':');
            if (separatorIndex > 0) {
                fileName = spectrumId.substring(0, separatorIndex);
            }
        }

        if (fileName == null || fileName.isEmpty()) {

            String[] fileNames = spectrumProvider.getOrderedFileNamesWithoutExtensions();

            if (fileNames != null && fileNames.length == 1) {
                fileName = fileNames[0];
            }
        }

        if (fileName == null || fileName.isEmpty()) {
            throw new IllegalArgumentException("Unable to infer the spectrum file name from the InstaNovo csv file.");
        }

        return IoUtil.removeExtension(fileName);
    }

    /**
     * Resolves the spectrum title.
     *
     * @param spectrumProvider the spectrum provider
     * @param spectrumFileName the spectrum file name without extension
     * @param spectrumId the spectrum id
     * @param scanNumber the scan number
     *
     * @return the spectrum title
     */
    private String getSpectrumTitle(SpectrumTitleLookup spectrumTitleLookup, String spectrumFileName, String spectrumId, String scanNumber) {

        String title = spectrumTitleLookup.getTitle(spectrumId);

        if (title != null) {
            return title;
        }

        if (spectrumId != null) {
            int separatorIndex = spectrumId.indexOf(':');

            if (separatorIndex >= 0 && separatorIndex < spectrumId.length() - 1) {

                title = spectrumTitleLookup.getTitle(spectrumId.substring(separatorIndex + 1));

                if (title != null) {
                    return title;
                }
            }
        }

        if (scanNumber != null && !scanNumber.isEmpty()) {

            title = spectrumTitleLookup.getTitle(scanNumber);

            if (title != null) {
                return title;
            }

            title = spectrumTitleLookup.getTitleForNumber(scanNumber);

            if (title != null) {
                return title;
            }
        }

        if (spectrumId != null) {
            int separatorIndex = spectrumId.indexOf(':');

            if (separatorIndex >= 0 && separatorIndex < spectrumId.length() - 1) {

                title = spectrumTitleLookup.getTitleAtIndex(spectrumId.substring(separatorIndex + 1));

                if (title != null) {
                    return title;
                }
            }
        }

        throw new IllegalArgumentException("Unable to match InstaNovo spectrum id '" + spectrumId + "' to a spectrum title in file '" + spectrumFileName + "'.");
    }

    /**
     * Returns the precursor charge.
     *
     * @param value the charge column value
     * @param lineNumber the line number
     * @param waitingHandler the waiting handler
     *
     * @return the charge, or null if the row should be skipped
     */
    private Integer getCharge(String value, int lineNumber, WaitingHandler waitingHandler) {

        String charge = value == null ? "" : value.trim();

        try {
            return Integer.parseInt(charge);
        } catch (NumberFormatException e) {

            if (waitingHandler != null) {
                waitingHandler.appendReport(
                        "Skipping InstaNovo csv line " + lineNumber + ": invalid precursor charge '" + charge + "'.",
                        true,
                        true
                );
            }

            return null;
        }
    }

    /**
     * Parses a peptide sequence with optional UniMod annotations.
     *
     * @param prediction the prediction
     * @param lineNumber the line number
     *
     * @return the parsed peptide
     */
    private ParsedPeptide parsePeptide(String prediction, int lineNumber) {

        StringBuilder sequence = new StringBuilder();
        ArrayList<ModificationMatch> modifications = new ArrayList<>();
        int lastResidueSite = 0;

        for (int i = 0; i < prediction.length(); i++) {

            char currentChar = prediction.charAt(i);

            if (currentChar == '[') {

                int endIndex = prediction.indexOf(']', i);

                if (endIndex < 0) {
                    throw new IllegalArgumentException("Invalid UniMod annotation in InstaNovo csv file at line " + lineNumber + ".");
                }

                String annotation = prediction.substring(i + 1, endIndex);
                Character previousResidue = lastResidueSite > 0 ? sequence.charAt(lastResidueSite - 1) : null;
                Character nextResidue = previousResidue == null ? getNextResidue(prediction, endIndex + 1) : null;
                UtilitiesModification modification = getUtilitiesModification(annotation, previousResidue, nextResidue, lastResidueSite);

                if (modification != null) {
                    modifications.add(new ModificationMatch(modification.name, modification.site));
                }

                i = endIndex;

            } else if (Character.isLetter(currentChar)) {

                sequence.append(Character.toUpperCase(currentChar));
                lastResidueSite = sequence.length();
            }
        }

        if (sequence.length() == 0) {
            throw new IllegalArgumentException("No peptide sequence found in InstaNovo csv file at line " + lineNumber + ".");
        }

        return new ParsedPeptide(sequence.toString(), modifications.toArray(new ModificationMatch[modifications.size()]));
    }

    /**
     * Maps InstaNovo UniMod annotations to Utilities modification names.
     *
     * @param annotation the annotation
     * @param previousResidue the preceding residue, null for N-terminal
     * annotations
     * @param nextResidue the next residue, null when unavailable
     * @param site the preceding residue site
     *
     * @return the Utilities modification, or null if unsupported
     */
    private UtilitiesModification getUtilitiesModification(String annotation, Character previousResidue, Character nextResidue, int site) {

        if (!annotation.toUpperCase().startsWith("UNIMOD:")) {
            return null;
        }

        String accession = annotation.substring("UNIMOD:".length());

        if ("1".equals(accession) && previousResidue == null) {
            return new UtilitiesModification("Acetylation of peptide N-term", 0);
        } else if ("4".equals(accession) && previousResidue != null && previousResidue == 'C') {
            return new UtilitiesModification("Carbamidomethylation of C", site);
        } else if ("5".equals(accession) && previousResidue == null) {
            return new UtilitiesModification("Carbamilation of protein N-term", 0);
        } else if ("7".equals(accession) && previousResidue != null) {
            if (previousResidue == 'N') {
                return new UtilitiesModification("Deamidation of N", site);
            } else if (previousResidue == 'Q') {
                return new UtilitiesModification("Deamidation of Q", site);
            } else if (previousResidue == 'R') {
                return new UtilitiesModification("Citrullination of R", site);
            }
        } else if ("35".equals(accession) && previousResidue != null) {
            if (previousResidue == 'M') {
                return new UtilitiesModification("Oxidation of M", site);
            } else if (previousResidue == 'P') {
                return new UtilitiesModification("Oxidation of P", site);
            } else if (previousResidue == 'K') {
                return new UtilitiesModification("Oxidation of K", site);
            } else if (previousResidue == 'C') {
                return new UtilitiesModification("Oxidation of C", site);
            } else if (previousResidue == 'N') {
                return new UtilitiesModification("Oxidation of N", site);
            }
        } else if ("21".equals(accession) && previousResidue != null) {
            if (previousResidue == 'S') {
                return new UtilitiesModification("Phosphorylation of S", site);
            } else if (previousResidue == 'T') {
                return new UtilitiesModification("Phosphorylation of T", site);
            } else if (previousResidue == 'Y') {
                return new UtilitiesModification("Phosphorylation of Y", site);
            }
        } else if ("385".equals(accession)) {
            if (previousResidue != null && previousResidue == 'N' && site > 0) {
                return new UtilitiesModification("Ammonia loss from N", site);
            } else if (previousResidue != null && previousResidue == 'C' && site == 1) {
                return new UtilitiesModification("Pyrolidone from carbamidomethylated C", site);
            } else if (previousResidue == null && nextResidue != null) {
                if (nextResidue == 'N') {
                    return new UtilitiesModification("Ammonia loss from N", 1);
                } else if (nextResidue == 'C') {
                    return new UtilitiesModification("Pyrolidone from carbamidomethylated C", 1);
                }
            }
        }

        return null;
    }

    /**
     * Returns the next residue in the prediction.
     *
     * @param prediction the prediction
     * @param startIndex the start index
     *
     * @return the next residue, or null
     */
    private Character getNextResidue(String prediction, int startIndex) {

        for (int i = startIndex; i < prediction.length(); i++) {

            char currentChar = prediction.charAt(i);

            if (Character.isLetter(currentChar)) {
                return Character.toUpperCase(currentChar);
            }
        }

        return null;
    }

    /**
     * Returns a value from a parsed CSV row.
     *
     * @param values the values
     * @param index the index
     *
     * @return the value
     */
    private String getValue(ArrayList<String> values, int index) {
        return index < values.size() ? values.get(index) : "";
    }

    /**
     * Returns indexes by lowercase header.
     *
     * @param headers the headers
     *
     * @return the indexes
     */
    private HashMap<String, Integer> getColumnIndexes(ArrayList<String> headers) {

        HashMap<String, Integer> result = new HashMap<>();

        for (int i = 0; i < headers.size(); i++) {
            result.put(headers.get(i).trim().toLowerCase(), i);
        }

        return result;
    }

    /**
     * Returns an optional column.
     *
     * @param columnIndexes the column indexes
     * @param columnNames the column names
     *
     * @return the column index, or -1
     */
    private int getOptionalColumn(HashMap<String, Integer> columnIndexes, String... columnNames) {

        for (String columnName : columnNames) {

            Integer columnIndex = columnIndexes.get(columnName.toLowerCase());

            if (columnIndex != null) {
                return columnIndex;
            }
        }

        return -1;
    }

    /**
     * Returns a required column.
     *
     * @param columnIndexes the column indexes
     * @param columnNames the column names
     *
     * @return the column index
     */
    private int getRequiredColumn(HashMap<String, Integer> columnIndexes, String... columnNames) {

        int columnIndex = getOptionalColumn(columnIndexes, columnNames);

        if (columnIndex < 0) {
            throw new IllegalArgumentException("Mandatory columns are missing in the InstaNovo csv file.");
        }

        return columnIndex;
    }

    /**
     * Parses a CSV line.
     *
     * @param line the line
     *
     * @return the values
     */
    private ArrayList<String> parseCsvLine(String line) {

        ArrayList<String> values = new ArrayList<>();
        StringBuilder currentValue = new StringBuilder();
        boolean inQuotes = false;

        for (int i = 0; i < line.length(); i++) {

            char currentChar = line.charAt(i);

            if (currentChar == '"') {
                if (inQuotes && i + 1 < line.length() && line.charAt(i + 1) == '"') {
                    currentValue.append('"');
                    i++;
                } else {
                    inQuotes = !inQuotes;
                }
            } else if (currentChar == ',' && !inQuotes) {
                values.add(currentValue.toString());
                currentValue.setLength(0);
            } else {
                currentValue.append(currentChar);
            }
        }

        values.add(currentValue.toString());

        return values;
    }

    /**
     * Spectrum title lookup cache for one spectrum file.
     */
    private static class SpectrumTitleLookup {

        /**
         * Titles indexed by exact and lower-case title.
         */
        private final HashMap<String, String> titles = new HashMap<>();
        /**
         * Titles indexed by scan or index number tokens parsed from the title.
         */
        private final HashMap<String, String> titleByNumber = new HashMap<>();
        /**
         * Titles in spectrum file order.
         */
        private final String[] orderedTitles;

        /**
         * Constructor.
         *
         * @param spectrumProvider the spectrum provider
         * @param spectrumFileName the spectrum file name without extension
         */
        private SpectrumTitleLookup(SpectrumProvider spectrumProvider, String spectrumFileName) {

            String[] spectrumTitles = spectrumProvider.getSpectrumTitles(spectrumFileName);

            if (spectrumTitles == null || spectrumTitles.length == 0) {
                throw new IllegalArgumentException("No spectra found for file '" + spectrumFileName + "'.");
            }

            orderedTitles = spectrumTitles;

            for (String title : spectrumTitles) {
                addTitle(title);
            }
        }

        /**
         * Adds a title.
         *
         * @param title the title
         */
        private void addTitle(String title) {

            if (title == null) {
                return;
            }

            titles.put(title, title);
            titles.put(title.toLowerCase(), title);

            Matcher matcher = TITLE_NUMBER_PATTERN.matcher(title);

            while (matcher.find()) {
                addNumber(matcher.group(1), title);
            }
        }

        /**
         * Adds a scan or index number.
         *
         * @param number the number
         * @param title the spectrum title
         */
        private void addNumber(String number, String title) {

            String normalizedNumber = normalizeNumber(number);

            if (normalizedNumber == null) {
                return;
            }

            if (titleByNumber.containsKey(normalizedNumber)
                    && !title.equals(titleByNumber.get(normalizedNumber))) {
                titleByNumber.put(normalizedNumber, null);
            } else {
                titleByNumber.put(normalizedNumber, title);
            }
        }

        /**
         * Returns a title matching the given title candidate.
         *
         * @param candidate the candidate
         *
         * @return the title, or null if not found
         */
        private String getTitle(String candidate) {

            if (candidate == null) {
                return null;
            }

            String trimmedCandidate = candidate.trim();

            if (trimmedCandidate.isEmpty()) {
                return null;
            }

            String result = titles.get(trimmedCandidate);

            if (result != null) {
                return result;
            }

            return titles.get(trimmedCandidate.toLowerCase());
        }

        /**
         * Returns a title matching the given scan or index number.
         *
         * @param candidate the candidate
         *
         * @return the title, or null if not found
         */
        private String getTitleForNumber(String candidate) {

            String normalizedNumber = normalizeNumber(candidate);

            return normalizedNumber == null ? null : titleByNumber.get(normalizedNumber);
        }

        /**
         * Returns a title by zero-based spectrum position.
         *
         * @param candidate the candidate index
         *
         * @return the title, or null if not found
         */
        private String getTitleAtIndex(String candidate) {

            String normalizedNumber = normalizeNumber(candidate);

            if (normalizedNumber == null) {
                return null;
            }

            int index;

            try {
                index = Integer.parseInt(normalizedNumber);
            } catch (NumberFormatException e) {
                return null;
            }

            return index >= 0 && index < orderedTitles.length ? orderedTitles[index] : null;
        }

        /**
         * Normalizes a positive integer string.
         *
         * @param number the number
         *
         * @return the normalized number
         */
        private String normalizeNumber(String number) {

            if (number == null) {
                return null;
            }

            String trimmedNumber = number.trim();

            if (trimmedNumber.isEmpty()) {
                return null;
            }

            for (int i = 0; i < trimmedNumber.length(); i++) {
                if (!Character.isDigit(trimmedNumber.charAt(i))) {
                    return null;
                }
            }

            int startIndex = 0;

            while (startIndex < trimmedNumber.length() - 1 && trimmedNumber.charAt(startIndex) == '0') {
                startIndex++;
            }

            return trimmedNumber.substring(startIndex);
        }
    }

    /**
     * Parsed peptide values.
     */
    private static class ParsedPeptide {

        /**
         * The bare sequence.
         */
        private final String sequence;
        /**
         * The variable modifications.
         */
        private final ModificationMatch[] modificationMatches;

        /**
         * Constructor.
         *
         * @param sequence the sequence
         * @param modificationMatches the modification matches
         */
        private ParsedPeptide(String sequence, ModificationMatch[] modificationMatches) {
            this.sequence = sequence;
            this.modificationMatches = modificationMatches;
        }
    }

    /**
     * Utilities modification mapping.
     */
    private static class UtilitiesModification {

        /**
         * The modification name.
         */
        private final String name;
        /**
         * The modification site.
         */
        private final int site;

        /**
         * Constructor.
         *
         * @param name the modification name
         * @param site the modification site
         */
        private UtilitiesModification(String name, int site) {
            this.name = name;
            this.site = site;
        }
    }
}
