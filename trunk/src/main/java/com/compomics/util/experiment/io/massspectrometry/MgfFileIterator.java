package com.compomics.util.experiment.io.massspectrometry;

import com.compomics.util.experiment.massspectrometry.MSnSpectrum;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 * An iterator of the spectra in an mgf file.
 *
 * @author Marc Vaudel
 */
public class MgfFileIterator {

    /**
     * The reader going through the file.
     */
    private BufferedReader br;
    /**
     * The next spectrum in the file.
     */
    private MSnSpectrum nextSpectrum = null;
    /**
     * The name of the mgf file.
     */
    private String mgfFileName;
    /**
     * The rank of the spectrum.
     */
    private int rank;

    /**
     * Constructor.
     *
     * @param mgfFile the file to go through
     *
     * @throws FileNotFoundException if a FileNotFoundException occurs
     * @throws IOException if an IOException occurs
     */
    public MgfFileIterator(File mgfFile) throws FileNotFoundException, IOException {
        mgfFileName = mgfFile.getName();
        br = new BufferedReader(new FileReader(mgfFile));
        nextSpectrum = MgfReader.getSpectrum(br, mgfFileName);
        rank = 1;
        if (nextSpectrum.getScanNumber() == null) {
            nextSpectrum.setScanNumber(rank + "");
        } else {
            while (nextSpectrum.getScanNumber().equals(++rank + ""));
        }
    }

    /**
     * Indicates whether the file contains another spectrum.
     *
     * @return a boolean indicating whether the file contains another spectrum
     */
    public boolean hasNext() {
        return nextSpectrum != null;
    }

    /**
     * Returns the next spectrum in the file.
     *
     * @return the next spectrum in the file
     *
     * @throws IOException if an IOException occurs
     */
    public MSnSpectrum next() throws IOException {

        MSnSpectrum currentSpectrum = nextSpectrum;
        nextSpectrum = MgfReader.getSpectrum(br, mgfFileName);

        if (nextSpectrum == null) {
            br.close();
        } else if (nextSpectrum.getScanNumber() == null) {
            nextSpectrum.setScanNumber(++rank + "");
        } else {
            while (nextSpectrum.getScanNumber().equals(++rank + ""));
        }

        return currentSpectrum;
    }
}
