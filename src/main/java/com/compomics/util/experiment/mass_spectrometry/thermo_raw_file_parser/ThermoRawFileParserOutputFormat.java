package com.compomics.util.experiment.mass_spectrometry.thermo_raw_file_parser;

/**
 * The mass spectrometry output formats supported by ThermoRawFileParser.
 *
 * @author Harald Barsnes
 */
public enum ThermoRawFileParserOutputFormat {

    /**
     * Mascot generic format.
     */
    mgf(0, "mgf", "Mascot generic format", ".mgf"),
    /**
     * mzML generic PSI format.
     */
    mzML(1, "mzML", "mzML generic PSI format", ".mzML"),
    /**
     * Indexed mzML generic PSI format.
     */
    mzMLIndexed(2, "mzML (indexed)", "mzML generic PSI format", ".mzML");

    /**
     * The index of the format.
     */
    public final int index;
    /**
     * The name of the format.
     */
    public final String name;
    /**
     * A brief description of the format.
     */
    public final String description;
    /**
     * The file name ending of the format.
     */
    public final String fileNameEnding;

    /**
     * Constructor.
     *
     * @param index the index of the format
     * @param name the name of the format
     * @param description the brief description of the format
     * @param fileNameEnding the file name ending of the format
     */
    private ThermoRawFileParserOutputFormat(int index, String name, String description, String fileNameEnding) {
        this.index = index;
        this.name = name;
        this.description = description;
        this.fileNameEnding = fileNameEnding;
    }

    @Override
    public String toString() {
        return name;
    }

    /**
     * Empty default constructor.
     */
    private ThermoRawFileParserOutputFormat() {
        index = 0;
        name = "";
        description = "";
        fileNameEnding = "";
    }
}
