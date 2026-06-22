package com.compomics.util.experiment.io.identification.idfilereaders;

import com.compomics.util.experiment.identification.Advocate;
import java.io.File;

/**
 * Reader for standalone InstaNovo+ CSV output.
 *
 * @author CompOmics
 */
public class InstaNovoPlusIdfileReader extends InstaNovoCsvIdfileReader {

    /**
     * The supported extension.
     */
    public static final String EXTENSION = ".instanovoplus.csv";

    /**
     * Default constructor for service loading.
     */
    public InstaNovoPlusIdfileReader() {
        this(null);
    }

    /**
     * Constructor.
     *
     * @param csvFile the CSV file
     */
    public InstaNovoPlusIdfileReader(File csvFile) {
        super(csvFile, Advocate.instanovoPlus, EXTENSION);
    }
}
