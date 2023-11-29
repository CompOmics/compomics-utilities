package com.compomics.util.experiment.biology.genes.ensembl;

import com.compomics.util.experiment.biology.taxonomy.mappings.EnsemblSpecies.EnsemblDivision;
import static com.compomics.util.experiment.biology.taxonomy.mappings.EnsemblSpecies.EnsemblDivision.fungi;
import static com.compomics.util.experiment.biology.taxonomy.mappings.EnsemblSpecies.EnsemblDivision.metazoa;
import static com.compomics.util.experiment.biology.taxonomy.mappings.EnsemblSpecies.EnsemblDivision.plants;
import static com.compomics.util.experiment.biology.taxonomy.mappings.EnsemblSpecies.EnsemblDivision.protists;

/**
 * Class for the handling of Ensembl versions.
 *
 * @author Marc Vaudel
 * @author Harald Barsnes
 */
public class EnsemblVersion {

    /**
     * Returns the current Ensembl version number.
     *
     * @param ensemblDivision the Ensembl genome division
     *
     * @return the current Ensembl version number
     */
    public static Integer getCurrentEnsemblVersion(EnsemblDivision ensemblDivision) {

        // @TODO: find a less hard coded way of finding the current ensembl versions!!!
        if (ensemblDivision != EnsemblDivision.vertebrates) {
            return 57;
        } else {
            return 110;
        }

        // the code below used to work but is not always updated when new ensembl versions are released
//        if (ensemblVersions == null) {
//            ensemblVersions = new HashMap<String, Integer>();
//        }
//        if (!ensemblVersions.containsKey(ensemblType)) {
//
//            try {
//                // get the current Ensembl version
//                URL url = new URL("http://www.biomart.org/biomart/martservice?type=registry");
//
//                BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
//
//                String inputLine;
//                boolean ensemblVersionFound = false;
//                String ensemblVersionAsText = "?";
//
//                while ((inputLine = in.readLine()) != null && !ensemblVersionFound) {
//                    if (inputLine.indexOf("database=\"" + ensemblType + "_mart_") != -1) {
//                        ensemblVersionAsText = inputLine.substring(inputLine.indexOf("database=\"" + ensemblType + "_mart_") + ("database=\"" + ensemblType + "_mart_").length());
//                        ensemblVersionAsText = ensemblVersionAsText.substring(0, ensemblVersionAsText.indexOf("\""));
//                        ensemblVersionFound = true;
//                    }
//                }
//
//                in.close();
//
//                if (ensemblVersionFound) {
//                    try {
//                        Integer ensemblVersion = new Integer(ensemblVersionAsText);
//                        ensemblVersions.put(ensemblType, ensemblVersion);
//                    } catch (NumberFormatException e) {
//                        e.printStackTrace();
//                    }
//                }
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//        }
//
//        return ensemblVersions.get(ensemblType);
    }

    /**
     * Returns the name of the Ensembl schema for BioMart queries.
     *
     * @param ensemblDivision the Ensembl division
     *
     * @return the name of the Ensembl schema for BioMart queries
     */
    public static String getEnsemblSchemaName(EnsemblDivision ensemblDivision) {

        if (ensemblDivision == null) {
            return "default";
        }

        switch (ensemblDivision) {

//            case fungi:
//                return "fungi_mart_" + getCurrentEnsemblVersion(ensemblGenomeDivision);
//            case plants:
//                return "plants_mart_" + getCurrentEnsemblVersion(ensemblGenomeDivision);
//            case protists:
//                return "protists_mart_" + getCurrentEnsemblVersion(ensemblGenomeDivision);
//            case metazoa:
//                return "metazoa_mart_" + getCurrentEnsemblVersion(ensemblGenomeDivision);
            case vertebrates:
                return "default";
            case fungi:
                return "fungi_mart";
            case plants:
                return "plants_mart";
            case protists:
                return "protists_mart";
            case metazoa:
                return "metazoa_mart";
            default:
                return "default";

        }

    }

}
