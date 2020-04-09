package com.compomics.cli.fasta;

import com.compomics.software.cli.CommandParameter;
import com.compomics.util.Util;
import com.compomics.util.experiment.io.biology.protein.FastaParameters;
import com.compomics.util.experiment.io.biology.protein.FastaSummary;
import com.compomics.util.io.IoUtil;
import java.io.File;
import java.io.IOException;
import java.util.Date;
import org.apache.commons.cli.CommandLine;

/**
 * This class gathers command line parameters for the parsing of FASTA files.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class FastaParametersInputBean {

    /**
     * The FASTA parsing parameters.
     */
    private FastaParameters fastaParameters;

    /**
     * Verifies the command line start parameters.
     *
     * @param aLine the command line to validate
     *
     * @return true if the startup was valid
     */
    public static boolean isValidStartup(CommandLine aLine) {

        if (aLine.getOptions().length == 0) {

            return false;

        }

        if (aLine.hasOption(FastaParametersCLIParams.SUFFIX.id)) {

            String arg = aLine.getOptionValue(FastaParametersCLIParams.SUFFIX.id);

            if (!CommandParameter.isInList(arg, arg, new String[]{"1", "2"})) {

                return false;

            }
        }

        return true;
    }

    /**
     * Parses all the arguments from a command line.
     *
     * @param aLine the command line
     * @param fastaFile the FASTA file to infer the parameters from if not
     * provided in the command line arguments
     *
     * @throws IOException if an error occurs while reading or writing a file.
     */
    public FastaParametersInputBean(CommandLine aLine, File fastaFile) throws IOException {

        FastaParameters tempFastaParameters = new FastaParameters();
        FastaParameters inferredParameters = null;
        
        if (aLine.hasOption(FastaParametersCLIParams.DECOY_FLAG.id)) {

            String arg = aLine.getOptionValue(FastaParametersCLIParams.DECOY_FLAG.id);

            if (!arg.equals("")) {
                tempFastaParameters.setDecoyFlag(arg);            
            }

        } else {

            inferredParameters = FastaParameters.inferParameters(fastaFile.getAbsolutePath());
            tempFastaParameters.setDecoyFlag(inferredParameters.getDecoyFlag());
            
        }
        
        if (aLine.hasOption(FastaParametersCLIParams.SUFFIX.id)) {

            String arg = aLine.getOptionValue(FastaParametersCLIParams.SUFFIX.id);

            if (arg.equals("1")) {

                tempFastaParameters.setDecoySuffix(false);

            } else {

                tempFastaParameters.setDecoySuffix(true);
            }
        } else {

            if (inferredParameters==null)
                inferredParameters = FastaParameters.inferParameters(fastaFile.getAbsolutePath());
            
            tempFastaParameters.setDecoySuffix(inferredParameters.isDecoySuffix());
        }

        FastaSummary fastaSummary = FastaSummary.getSummary(fastaFile.getAbsolutePath(), tempFastaParameters, null);

        if (aLine.hasOption(FastaParametersCLIParams.NAME.id)) {

            String arg = aLine.getOptionValue(FastaParametersCLIParams.NAME.id);
            fastaSummary.setName(arg);

        } else {

            String fileName = IoUtil.removeExtension(fastaFile.getName());
            fastaSummary.setName(fileName);

        }

        if (aLine.hasOption(FastaParametersCLIParams.DESCRIPTION.id)) {

            String arg = aLine.getOptionValue(FastaParametersCLIParams.DESCRIPTION.id);
            fastaSummary.setDescription(arg);

        } else {

            String fileName = IoUtil.removeExtension(fastaFile.getName());
            fastaSummary.setDescription(fileName);

        }

        if (aLine.hasOption(FastaParametersCLIParams.VERSION.id)) {

            String arg = aLine.getOptionValue(FastaParametersCLIParams.VERSION.id);
            fastaSummary.setVersion(arg);

        } else {

            String fileVersion = new Date(fastaFile.lastModified()).toString();
            fastaSummary.setName(fileVersion);

        }

        FastaSummary.saveSummary(fastaFile.getAbsolutePath(), fastaSummary);

        this.fastaParameters = tempFastaParameters;
    }

    /**
     * Returns the FASTA parameters as parsed from the command line.
     *
     * @return the FASTA parameters as parsed from the command line
     */
    public FastaParameters getFastaParameters() {
        return fastaParameters;
    }
}
