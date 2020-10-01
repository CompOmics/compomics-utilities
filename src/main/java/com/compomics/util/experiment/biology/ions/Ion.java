package com.compomics.util.experiment.biology.ions;

import com.compomics.util.experiment.biology.atoms.AtomChain;
import com.compomics.util.experiment.biology.ions.impl.PrecursorIon;
import com.compomics.util.experiment.biology.ions.impl.Glycan;
import com.compomics.util.experiment.biology.ions.impl.PeptideFragmentIon;
import com.compomics.util.experiment.biology.ions.impl.TagFragmentIon;
import com.compomics.util.experiment.biology.ions.impl.ReporterIon;
import com.compomics.util.experiment.biology.ions.impl.ElementaryIon;
import com.compomics.util.experiment.biology.ions.impl.RelatedIon;
import com.compomics.util.experiment.biology.ions.impl.ImmoniumIon;
import com.compomics.util.experiment.biology.aminoacids.AminoAcid;
import com.compomics.util.experiment.personalization.ExperimentObject;
import com.compomics.util.pride.CvTerm;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.stream.Collectors;

/**
 * This class models an ion.
 *
 * @author Marc Vaudel
 */
public abstract class Ion extends ExperimentObject {

    /**
     * Empty default constructor
     */
    public Ion() {
    }

    /**
     * Serial number for backward compatibility.
     */
    static final long serialVersionUID = -1505719074403886934L;
    /**
     * Cache for the neutral losses as string.
     */
    private String neutralLossesAsString = null;

    /**
     * An enumerator of the supported ion types.
     */
    public enum IonType {

        /**
         * Identifier for a peptide fragment ion.
         */
        PEPTIDE_FRAGMENT_ION(0),
        /**
         * A tag fragment ion
         */
        TAG_FRAGMENT_ION(1),
        /**
         * Identifier for an MH ion. The number of H is not represented here.
         */
        PRECURSOR_ION(2),
        /**
         * Identifier for an immonium ion.
         */
        IMMONIUM_ION(3),
        /**
         * Identifier for a reporter ion.
         */
        REPORTER_ION(4),
        /**
         * Identifier for a glycan.
         */
        GLYCAN(5),
        /**
         * Identifier for an elementary ion.
         */
        ELEMENTARY_ION(6),
        /**
         * Identifier for an unknown ion.
         */
        UNKNOWN(7),
        /**
         * Identifier for a related ion.
         */
        RELATED_ION(8);

        /**
         * The index of the type.
         */
        public final int index;

        /**
         * Constructor.
         *
         * @param index the index of the type
         */
        private IonType(int index) {
            this.index = index;
        }
    }
    /**
     * Type of ion.
     */
    protected IonType type = IonType.UNKNOWN;
    /**
     * The theoretic mass.
     * 
     * @deprecated use the double value instead.
     */
    protected Double theoreticMass;
    /**
     * The theoretic mass.
     */
    protected double theoreticMass1;
    /**
     * The atomic composition of the ion.
     */
    protected AtomChain atomChain;

    /**
     * Returns the name of the ion. The name should be short enough to be
     * displayed on a spectrum.
     *
     * @return the name of the ion
     */
    public abstract String getName();

    /**
     * Returns the CV term adapted to the fragment ion. Null if none
     * corresponding.
     *
     * @return the CV term adapted to the fragment ion. Null if none
     * corresponding
     */
    public abstract CvTerm getPrideCvTerm();

    /**
     * Returns the CV term adapted to the fragment ion. Null if none
     * corresponding.
     *
     * @return the CV term adapted to the fragment ion. Null if none
     * corresponding
     */
    public abstract CvTerm getPsiMsCvTerm();

    /**
     * Returns the ion subtype.
     *
     * @return the ion subtype as integer
     */
    public abstract int getSubType();

    /**
     * Returns the subtype as string.
     *
     * @return the subtype as string
     */
    public abstract String getSubTypeAsString();

    /**
     * Returns an array of possible subtypes.
     *
     * @param ionType an array of possible subtypes
     * 
     * @return an array of possible subtypes
     */
    public static int[] getPossibleSubtypes(IonType ionType) {
        switch (ionType) {
            case ELEMENTARY_ION:
                return ElementaryIon.getPossibleSubtypes();
            case GLYCAN:
                return Glycan.getPossibleSubtypes();
            case IMMONIUM_ION:
                return ImmoniumIon.getPossibleSubtypes();
            case PEPTIDE_FRAGMENT_ION:
                return PeptideFragmentIon.getPossibleSubtypes();
            case TAG_FRAGMENT_ION:
                return TagFragmentIon.getPossibleSubtypes();
            case PRECURSOR_ION:
                return PrecursorIon.getPossibleSubtypes();
            case REPORTER_ION:
                return ReporterIon.getPossibleSubtypes();
            case RELATED_ION:
                return RelatedIon.getPossibleSubtypes();
            default:
                throw new UnsupportedOperationException("Not supported yet.");
        }
    }

    /**
     * Returns a hashset of possible subtypes.
     *
     * @param ionType a hashset of possible subtypes
     * 
     * @return a hashset of possible subtypes
     */
    public static HashSet<Integer> getPossibleSubtypesAsSet(IonType ionType) {
        
        int[] possibleSubtypes = getPossibleSubtypes(ionType);
        
        return Arrays.stream(possibleSubtypes)
                .boxed()
                .collect(Collectors.toCollection(HashSet::new));
        
    }

    /**
     * Returns the possible neutral losses of this ion type. An empty list if
     * none.
     *
     * @return the possible neutral losses of this ion type
     */
    public abstract NeutralLoss[] getNeutralLosses();

    /**
     * Indicates whether the ion has a neutral loss.
     *
     * @return a boolean indicating whether the ion has a neutral loss
     */
    public boolean hasNeutralLosses() {
        switch (type) {
            case PEPTIDE_FRAGMENT_ION:
            case TAG_FRAGMENT_ION:
            case PRECURSOR_ION:
                NeutralLoss[] neutralLosses = getNeutralLosses();
                return neutralLosses != null && neutralLosses.length > 0;
            default:
                return false;
        }
    }

    /**
     * Returns a boolean indicating whether the ion is the same as another ion.
     *
     * @param anotherIon the other ion
     * @return a boolean indicating whether the ion is the same as another ion
     */
    public abstract boolean isSameAs(Ion anotherIon);

    /**
     * Returns the neutral loss (if any), the empty string if no loss.
     *
     * @return the neutral loss
     */
    public String getNeutralLossesAsString() {
        if (neutralLossesAsString == null) {
            neutralLossesAsString = getNeutralLossesAsString(getNeutralLosses());
        }
        return neutralLossesAsString;
    }

    /**
     * Returns the neutral loss (if any), the empty string if no loss.
     *
     * @param neutralLosses the neutral loss (if any)
     * @return the neutral loss
     */
    public static String getNeutralLossesAsString(NeutralLoss[] neutralLosses) {
        if (neutralLosses == null) {
            return "";
        }
        ArrayList<String> names = new ArrayList<>(neutralLosses.length);
        for (NeutralLoss neutralLoss : neutralLosses) {
            names.add(neutralLoss.name);
        }
        Collections.sort(names);
        StringBuilder result = new StringBuilder(4 * neutralLosses.length);
        for (String name : names) {
            result.append('-').append(name);
        }
        return result.toString();
    }

    /**
     * Returns the theoretic mass, from the atomic composition if available,
     * from the theoreticMass field otherwise.
     *
     * @return the theoretic mass
     */
    public double getTheoreticMass() {
        if (atomChain != null) {
            return atomChain.getMass();
        }
        return theoreticMass1;
    }

    /**
     * Returns the m/z expected for this ion at the given charge.
     *
     * @param charge the charge of interest
     *
     * @return the m/z expected for this ion
     */
    public double getTheoreticMz(Integer charge) {
        double protonMass = ElementaryIon.proton.getTheoreticMass();
        double mz = getTheoreticMass() + protonMass;
        if (charge > 1) {
            mz = (mz + (charge - 1) * protonMass) / charge;
        }
        return mz;
    }

    /**
     * Returns the atomic composition.
     *
     * @return the atomic composition
     */
    public AtomChain getAtomicComposition() {
        return atomChain;
    }

    /**
     * Returns the atomic composition.
     *
     * @param atomChain the atomic composition
     */
    public void setAtomicComposition(AtomChain atomChain) {
        this.atomChain = atomChain;
    }

    /**
     * Sets a new theoretic mass.
     *
     * @param theoreticMass a new theoretic mass
     */
    public void setTheoreticMass(double theoreticMass) {
        this.theoreticMass1 = theoreticMass;
    }

    /**
     * Returns the ion type.
     *
     * @return the ion type
     */
    public IonType getType() {
        return type;
    }

    /**
     * Returns the implemented ion types.
     *
     * @return the implemented ion types
     */
    public static ArrayList<IonType> getImplementedIonTypes() {
        ArrayList<IonType> result = new ArrayList<>();
        result.add(IonType.ELEMENTARY_ION);
        result.add(IonType.GLYCAN);
        result.add(IonType.IMMONIUM_ION);
        result.add(IonType.PEPTIDE_FRAGMENT_ION);
        result.add(IonType.TAG_FRAGMENT_ION);
        result.add(IonType.PRECURSOR_ION);
        result.add(IonType.REPORTER_ION);
        result.add(IonType.RELATED_ION);
        return result;
    }

    /**
     * Returns the type of ion as string.
     *
     * @return the type of ion as string
     */
    public String getTypeAsString() {
        return getTypeAsString(type);
    }

    /**
     * Returns the type of ion as string.
     *
     * @param type the type of ion as string
     * @return the type of ion as string
     */
    public static String getTypeAsString(IonType type) {
        switch (type) {
            case PEPTIDE_FRAGMENT_ION:
                return "Peptide fragment ion";
            case TAG_FRAGMENT_ION:
                return "Tag fragment ion";
            case PRECURSOR_ION:
                return "Precursor ion";
            case IMMONIUM_ION:
                return "Immonium ion";
            case REPORTER_ION:
                return "Reporter ion";
            case GLYCAN:
                return "Glycan";
            case ELEMENTARY_ION:
                return "Elementary ion";
            case RELATED_ION:
                return "Related ion";
            case UNKNOWN:
                return "Unknown ion type";
            default:
                throw new UnsupportedOperationException("No name for ion type " + type + ".");
        }
    }

    /**
     * Convenience method returning a generic ion based on the given ion type.
     *
     * @param ionType the ion type
     * @param subType the ion subtype
     * @param neutralLosses the neutral losses. Null list if none.
     * @return a generic ion
     */
    public static Ion getGenericIon(IonType ionType, int subType, NeutralLoss[] neutralLosses) {
        switch (ionType) {
            case ELEMENTARY_ION:
                return new ElementaryIon("new ElementaryIon", 0.0, subType);
            case GLYCAN:
                return new Glycan("new Glycan", "new Glycan");
            case IMMONIUM_ION:
                return ImmoniumIon.getImmoniumIon(subType);
            case PEPTIDE_FRAGMENT_ION:
                return new PeptideFragmentIon(subType, neutralLosses);
            case TAG_FRAGMENT_ION:
                return new TagFragmentIon(subType, neutralLosses);
            case PRECURSOR_ION:
                return new PrecursorIon(neutralLosses);
            case REPORTER_ION:
                return ReporterIon.getReporterIon(subType);
            case RELATED_ION:
                return new RelatedIon(AminoAcid.A, AtomChain.getAtomChain("H"), -1, false);
            default:
                throw new UnsupportedOperationException("No generic constructor for " + getTypeAsString(ionType) + ".");
        }
    }

    /**
     * Convenience method returning a generic ion based on the given ion type
     * without neutral losses.
     *
     * @param ionType the ion type
     * @param subType the ion subtype
     * @return a generic ion
     */
    public static Ion getGenericIon(IonType ionType, int subType) {
        switch (ionType) {
            case ELEMENTARY_ION:
                return new ElementaryIon("new ElementaryIon", 0.0, subType);
            case GLYCAN:
                return new Glycan("new Glycon", "new Glycon");
            case IMMONIUM_ION:
                return ImmoniumIon.getImmoniumIon(subType);
            case PEPTIDE_FRAGMENT_ION:
                return new PeptideFragmentIon(subType);
            case TAG_FRAGMENT_ION:
                return new TagFragmentIon(subType);
            case PRECURSOR_ION:
                return new PrecursorIon();
            case REPORTER_ION:
                return ReporterIon.getReporterIon(subType);
            case RELATED_ION:
                return new RelatedIon(AminoAcid.A, AtomChain.getAtomChain("H"), -1, false);
            default:
                throw new UnsupportedOperationException("No generic constructor for " + getTypeAsString(ionType) + ".");
        }
    }
}
