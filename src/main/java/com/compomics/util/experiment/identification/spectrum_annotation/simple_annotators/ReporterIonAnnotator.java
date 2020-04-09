package com.compomics.util.experiment.identification.spectrum_annotation.simple_annotators;

import com.compomics.util.experiment.biology.ions.impl.ElementaryIon;
import com.compomics.util.experiment.biology.ions.impl.ReporterIon;
import com.compomics.util.experiment.identification.matches.IonMatch;
import com.compomics.util.experiment.mass_spectrometry.spectra.Peak;
import com.compomics.util.experiment.mass_spectrometry.indexes.SpectrumIndex;
import java.util.ArrayList;

/**
 * Annotator for reporter ions.
 *
 * @author Marc Vaudel
 */
public class ReporterIonAnnotator {

    /**
     * Empty default constructor
     */
    public ReporterIonAnnotator() {
        reporterIonsMz = null;
        reporterIons = null;
    }

    /**
     * Array of the m/z of the reporter ions to annotate.
     */
    private final double[] reporterIonsMz;

    /**
     * Array of the reporter ions to annotate.
     */
    private final ReporterIon[] reporterIons;

    /**
     * Constructor.
     *
     * @param reporterIons array of the reporter ions to annotate
     */
    public ReporterIonAnnotator(
            ReporterIon[] reporterIons
    ) {

        this.reporterIons = reporterIons;
        this.reporterIonsMz = new double[reporterIons.length];

        for (int i = 0; i < reporterIons.length; i++) {

            reporterIonsMz[i] = reporterIons[i].getTheoreticMass() + ElementaryIon.proton.getTheoreticMass();

        }
    }

    /**
     * Returns the ions matched in the given spectrum.
     *
     * @param spectrumIndex the index of the spectrum
     *
     * @return the ions matched in the given spectrum
     */
    public ArrayList<IonMatch> getIonMatches(
            SpectrumIndex spectrumIndex
    ) {

        ArrayList<IonMatch> results = new ArrayList<>(reporterIons.length);

        for (int i = 0; i < reporterIons.length; i++) {

            double ionMz = reporterIonsMz[i];
            int[] indexes = spectrumIndex.getMatchingPeaks(ionMz);

            if (indexes.length > 0) {

                ReporterIon ion = reporterIons[i];

                for (int index : indexes) {

                    results.add(
                            new IonMatch(
                                    spectrumIndex.mzArray[index],
                                    spectrumIndex.intensityArray[index],
                                    ion,
                                    1
                            )
                    );
                }
            }
            }

            return results;
            
        }
    }
