package com.compomics.util.experiment.io.identification.idfilereaders;

import com.compomics.util.experiment.identification.Advocate;
import java.io.File;

/**
 * Reader for InstaNovo transformer-only CSV output.
 *
 * @author CompOmics
 */
public class InstaNovoIdfileReader extends InstaNovoCsvIdfileReader {

    /**
     * The supported extension.
     */
    public static final String EXTENSION = ".instanovo.csv";

    /**
     * Default constructor for service loading.
     */
    public InstaNovoIdfileReader() {
        this(null);
    }

    /**
     * Constructor.
     *
     * @param csvFile the CSV file
     */
    public InstaNovoIdfileReader(File csvFile) {
        super(csvFile, Advocate.instanovo, EXTENSION);
    }
}
