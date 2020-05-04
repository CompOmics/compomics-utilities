package com.compomics.util.experiment.mass_spectrometry;

import com.compomics.util.experiment.mass_spectrometry.spectra.Precursor;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import java.util.HashMap;

/**
 * Interface for objects providing spectra.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public interface SpectrumProvider extends AutoCloseable {

    /**
     * Returns the spectrum with the given title in the given file.
     *
     * @param fileName The name of the spectrum file without file extension.
     * @param spectrumTitle The title of the spectrum.
     *
     * @return The spectrum with the given title in the given file.
     */
    public Spectrum getSpectrum(
            String fileName,
            String spectrumTitle
    );

    /**
     * Returns the precursor. Null if none.
     *
     * @param fileName The name of the spectrum file without file extension.
     * @param spectrumTitle The title of the spectrum.
     *
     * @return The precursor.
     */
    public Precursor getPrecursor(
            String fileName,
            String spectrumTitle
    );

    /**
     * Returns the measured precursor m/z. NaN if none.
     *
     * @param fileName The name of the spectrum file without file extension.
     * @param spectrumTitle The title of the spectrum.
     *
     * @return The measured precursor m/z.
     */
    public double getPrecursorMz(
            String fileName,
            String spectrumTitle
    );

    /**
     * Returns the precursor RT window. NaN if none.
     *
     * @param fileName The name of the spectrum file without file extension.
     * @param spectrumTitle The title of the spectrum.
     *
     * @return The precursor RT.
     */
    public double getPrecursorRt(
            String fileName,
            String spectrumTitle
    );

    /**
     * Returns the spectrum peaks.
     *
     * @param fileName The name of the spectrum file without file extension.
     * @param spectrumTitle The title of the spectrum.
     *
     * @return The peaks.
     */
    public double[][] getPeaks(
            String fileName,
            String spectrumTitle
    );

    /**
     * Returns the minimum precursor m/z in a given file.
     *
     * @param fileName The name of the spectrum file without file extension.
     *
     * @return The minimum precursor m/z in a given file.
     */
    public double getMinPrecMz(
            String fileName
    );

    /**
     * Returns the maximum precursor m/z in a given file.
     *
     * @param fileName The name of the spectrum file without file extension.
     *
     * @return The maximum precursor m/z in a given file.
     */
    public double getMaxPrecMz(
            String fileName
    );

    /**
     * Returns the maximum precursor intensity in a given file.
     *
     * @param fileName The name of the spectrum file without file extension.
     *
     * @return The maximum precursor intensity in a given file.
     */
    public double getMaxPrecInt(
            String fileName
    );

    /**
     * Returns the maximum precursor RT in a given file.
     *
     * @param fileName The name of the spectrum file without file extension.
     *
     * @return The maximum precursor RT in a given file.
     */
    public double getMaxPrecRT(
            String fileName
    );

    /**
     * Returns the minimum precursor m/z among all files.
     *
     * @return The minimum precursor m/z among all files.
     */
    public double getMinPrecMz();

    /**
     * Returns the maximum precursor m/z among all files.
     *
     * @return The maximum precursor m/z among all files.
     */
    public double getMaxPrecMz();

    /**
     * Returns the maximum precursor intensity among all files.
     *
     * @return The maximum precursor intensity among all files.
     */
    public double getMaxPrecInt();

    /**
     * Returns the maximum precursor RT among all files.
     *
     * @return The maximum precursor RT among all files.
     */
    public double getMaxPrecRT();

    /**
     * Returns the spectrum file names without file extensions.
     *
     * @return The spectrum file names without file extensions.
     */
    public String[] getFileNames();

    /**
     * Returns the spectrum titles for the given mass spectrometry file name.
     *
     * @param fileName The mass spectrometry file name without file extension.
     *
     * @return The spectrum titles as array.
     */
    public String[] getSpectrumTitles(
            String fileName
    );

    /**
     * Returns the absolute path to the original mass spec file containing the
     * spectra in a map indexed by file name without file extension.
     *
     * @return The absolute path to the original mass spec file containing the
     * spectra in a map indexed by file name.
     */
    public HashMap<String, String> getFilePaths();

    /**
     * Returns the absolute path to the cms file indexed by ms file name without
     * file extension. Null if none.
     *
     * @return The absolute path to the cms file indexed by ms file name.
     */
    public HashMap<String, String> getCmsFilePaths();

    @Override
    public void close();

}
