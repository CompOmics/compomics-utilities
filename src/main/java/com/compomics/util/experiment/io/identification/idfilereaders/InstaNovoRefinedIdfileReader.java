package com.compomics.util.experiment.io.identification.idfilereaders;

import com.compomics.util.experiment.identification.Advocate;
import java.io.File;

/**
 * Reader for InstaNovo predictions refined by InstaNovo+.
 *
 * @author CompOmics
 */
public class InstaNovoRefinedIdfileReader extends InstaNovoCsvIdfileReader {

    /**
     * The supported extension.
     */
    public static final String EXTENSION = ".instanovo.refined.csv";

    /**
     * Default constructor for service loading.
     */
    public InstaNovoRefinedIdfileReader() {
        this(null);
    }

    /**
     * Constructor.
     *
     * @param csvFile the CSV file
     */
    public InstaNovoRefinedIdfileReader(File csvFile) {
        super(csvFile, Advocate.instanovoRefined, EXTENSION);
    }
}
