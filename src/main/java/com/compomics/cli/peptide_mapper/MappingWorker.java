package com.compomics.cli.peptide_mapper;

import com.compomics.util.experiment.biology.aminoacids.sequence.AminoAcidSequence;
import com.compomics.util.experiment.identification.protein_inference.FastaMapper;
import com.compomics.util.experiment.identification.protein_inference.PeptideProteinMapping;
import com.compomics.util.gui.waiting.waitinghandlers.WaitingHandlerCLIImpl;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import com.compomics.util.experiment.identification.amino_acid_tags.Tag;
import com.compomics.util.experiment.identification.protein_inference.fm_index.FMIndex;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;
import com.compomics.util.experiment.identification.utils.PeptideUtils;
import com.compomics.util.io.flat.SimpleFileReader;
import com.compomics.util.io.flat.SimpleFileWriter;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.threading.SimpleSemaphore;

/**
 * MappingWorker.
 *
 * @author dominik.kopczynski
 */
public class MappingWorker implements Runnable {

    private final WaitingHandlerCLIImpl waitingHandlerCLIImpl;
    private final FastaMapper peptideMapper;
    private final SequenceMatchingParameters sequenceMatchingPreferences;
    private final SimpleSemaphore bufferMutex;
    private final SimpleFileReader reader;
    private final SimpleFileWriter writer;
    private static final int NUM_READS = 1000;
    private final boolean flanking;
    private final boolean peptideMapping;
    public Exception exception = null;

    public MappingWorker(WaitingHandlerCLIImpl waitingHandlerCLIImpl,
            FastaMapper peptideMapper,
            IdentificationParameters identificationParameters,
            SimpleFileReader reader,
            SimpleSemaphore bufferMutex,
            SimpleFileWriter writer,
            boolean peptideMapping
    ) {
        this.waitingHandlerCLIImpl = waitingHandlerCLIImpl;
        this.peptideMapper = peptideMapper;
        this.sequenceMatchingPreferences = identificationParameters.getSequenceMatchingParameters();
        this.reader = reader;
        this.bufferMutex = bufferMutex;
        this.writer = writer;
        this.flanking = identificationParameters.getSearchParameters().getFlanking();
        this.peptideMapping = peptideMapping;
    }

    public String flanking(PeptideProteinMapping peptideProteinMapping, FastaMapper peptideMapper) {

        String peptide = peptideProteinMapping.getPeptideSequence();
        String accession = peptideProteinMapping.getProteinAccession();
        int peptideLength = peptide.length();

        char prefixChar = ((FMIndex) peptideMapper).prefixCharacter(accession, peptideProteinMapping.fmIndexPosition);
        if (prefixChar != FMIndex.DELIMITER) {
            peptide = Character.toString(prefixChar) + "." + peptide;
        } else {
            peptide = "-" + peptide;
        }

        char suffixChar = ((FMIndex) peptideMapper).suffixCharacter(accession, peptideProteinMapping.fmIndexPosition, peptideLength + 1);
        if (suffixChar != FMIndex.DELIMITER) {
            peptide += "." + Character.toString(suffixChar);

        } else {
            peptide += "-";
        }

        return peptide;
    }

    @Override
    public void run() {

        ArrayList<String> rows = new ArrayList<>();
        ArrayList<String> outputData = new ArrayList<>();

        while (true) {
            rows.clear();
            outputData.clear();

            // read input file batch wise
            try {
                String row;

                bufferMutex.acquire();

                int i = 0;

                while (!waitingHandlerCLIImpl.isRunCanceled() && i++ < NUM_READS && (row = reader.readLine()) != null) {

                    row = row.trim();

                    if (!row.isEmpty()) {

                        rows.add(row);

                    }
                }

                bufferMutex.release();

                if (waitingHandlerCLIImpl.isRunCanceled() || rows.isEmpty()) {
                    break;
                }

            } catch (Exception e) {
                waitingHandlerCLIImpl.setRunCanceled();
                exception = new IOException("Error: cound not open input list.\n\n" + e);
                return;
            }

            // map peptides sequences
            if (peptideMapping) {
                for (String inputPeptide : rows) {
                    if (waitingHandlerCLIImpl.isRunCanceled()) {
                        break;
                    }
                    for (char c : inputPeptide.toCharArray()) {
                        if (!(((int) 'A' <= c && c <= (int) 'Z') || ((int) 'a' <= c && c <= (int) 'z'))) {
                            waitingHandlerCLIImpl.setRunCanceled();
                            exception = new RuntimeException("Error: invalid character in line '" + inputPeptide + "' -> '" + (char) c + "'.");
                            return;
                        }
                    }

                    try {
                        for (PeptideProteinMapping peptideProteinMapping : peptideMapper.getProteinMapping(inputPeptide.toUpperCase(), sequenceMatchingPreferences)) {
                            String peptide = peptideProteinMapping.getPeptideSequence();

                            String accession = peptideProteinMapping.getProteinAccession();
                            int startIndex = peptideProteinMapping.getIndex() + 1;
                            if (flanking) {
                                peptide = flanking(peptideProteinMapping, peptideMapper);
                            }

                            String modifications = "";

                            if (peptideProteinMapping.getVariableModifications() != null) {
                                modifications = "," + PeptideUtils.getVariableModificationsAsString(peptideProteinMapping.getVariableModifications());
                            }

                            outputData.add(
                                    String.join(",", peptide, accession, startIndex + modifications)
                            );
                        }
                        waitingHandlerCLIImpl.increaseSecondaryProgressCounter();
                    } catch (Exception e) {
                        exception = new RuntimeException("An error occurred during the mapping of '" + inputPeptide + "'\n\n" + e);
                        waitingHandlerCLIImpl.setRunCanceled();
                    }
                }

            } else {
                for (String tagString : rows) {
                    if (waitingHandlerCLIImpl.isRunCanceled()) {
                        break;
                    }

                    Tag tag = new Tag();
                    for (String part : tagString.split(",")) {

                        if (Pattern.matches("[a-zA-Z]+", part)) {
                            tag.addAminoAcidSequence(new AminoAcidSequence(part));
                        } else {
                            try {
                                double mass = Double.parseDouble(part);
                                tag.addMassGap(mass);
                            } catch (NumberFormatException e) {
                                waitingHandlerCLIImpl.setRunCanceled();
                                exception = new RuntimeException("Error: line contains no valid tag: '" + tagString + "'.\n\n" + e);
                                return;
                            }
                        }
                    }

                    try {
                        for (PeptideProteinMapping peptideProteinMapping : peptideMapper.getProteinMapping(tag, sequenceMatchingPreferences)) {
                            String peptide = peptideProteinMapping.getPeptideSequence();

                            String accession = peptideProteinMapping.getProteinAccession();
                            int startIndex = peptideProteinMapping.getIndex() + 1;
                            if (flanking) {
                                peptide = flanking(peptideProteinMapping, peptideMapper);
                            }

                            String modifications = "";

                            if (peptideProteinMapping.getVariableModifications() != null) {
                                modifications = "," + PeptideUtils.getVariableModificationsAsString(peptideProteinMapping.getVariableModifications());
                            }

                            outputData.add(
                                    String.join(",", tagString, peptide, accession, startIndex + modifications)
                            );
                        }
                        waitingHandlerCLIImpl.increaseSecondaryProgressCounter();
                    } catch (Exception e) {
                        exception = new RuntimeException("An error occurred during the mapping of '" + tagString + "'\n\n" + e);
                        waitingHandlerCLIImpl.setRunCanceled();
                    }
                }
            }

            // write out processed batch
            try {
                for (String output : outputData) {
                    if (waitingHandlerCLIImpl.isRunCanceled()) {
                        break;
                    }
                    writer.writeLine(output);
                }
            } catch (Exception e) {
                exception = new IOException("Error: could not write into file.\n\n" + e);
                waitingHandlerCLIImpl.setRunCanceled();
                return;
            }
        }
    }
}
