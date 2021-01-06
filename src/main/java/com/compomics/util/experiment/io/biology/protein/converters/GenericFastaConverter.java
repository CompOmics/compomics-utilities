package com.compomics.util.experiment.io.biology.protein.converters;

import com.compomics.util.experiment.biology.proteins.Protein;
import com.compomics.util.experiment.io.biology.protein.Header;
import com.compomics.util.experiment.io.biology.protein.iterators.FastaIterator;
import com.compomics.util.waiting.WaitingHandler;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * This converter writes a FASTA file with standardized headers.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class GenericFastaConverter {

    /**
     * Appends decoy sequences to the provided FASTA file.
     *
     * @param fastaIn the FASTA file to read
     * @param fastaOut the FASTA file to write
     * @param waitingHandler a handler to allow canceling the import and
     * displaying progress
     *
     * @throws IOException exception thrown whenever an error happened while
     * reading or writing a FASTA file
     */
    public static void convertFile(File fastaIn, File fastaOut, WaitingHandler waitingHandler) throws IOException {

        if (!fastaOut.getParentFile().exists()) {
            fastaOut.getParentFile().mkdirs();
        }
        
        FastaIterator fastaIterator = new FastaIterator(fastaIn);

        try (BufferedWriter bw = new BufferedWriter(new FileWriter(fastaOut))) {

            Protein protein;
            while ((protein = fastaIterator.getNextProtein()) != null) {

                Header header = fastaIterator.getLastHeader();
                String genericHeader = header.asGenericHeader();

                bw.write(genericHeader);
                bw.newLine();
                bw.write(protein.getSequence());
                bw.newLine();
                bw.newLine();

                if (waitingHandler != null) {

                    if (waitingHandler.isRunCanceled()) {

                        return;

                    }

                    waitingHandler.increaseSecondaryProgressCounter();

                }
            }
        }
    }
}
