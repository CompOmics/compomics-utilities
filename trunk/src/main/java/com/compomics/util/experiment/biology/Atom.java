package com.compomics.util.experiment.biology;

import com.compomics.util.experiment.biology.atoms.Carbon;
import com.compomics.util.experiment.biology.atoms.Helium;
import com.compomics.util.experiment.biology.atoms.Hydrogen;
import com.compomics.util.experiment.biology.atoms.Iodine;
import com.compomics.util.experiment.biology.atoms.Lithium;
import com.compomics.util.experiment.biology.atoms.Nitrogen;
import com.compomics.util.experiment.biology.atoms.Oxygen;
import com.compomics.util.experiment.biology.atoms.Phosphorus;
import com.compomics.util.experiment.biology.atoms.Selenium;
import com.compomics.util.experiment.biology.atoms.Sodium;
import com.compomics.util.experiment.biology.atoms.Sulfur;
import com.compomics.util.experiment.personalization.ExperimentObject;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * This interface contains information about atoms.
 *
 * @author Marc Vaudel
 */
public abstract class Atom extends ExperimentObject {

    /**
     * The version UID for Serialization/Deserialization compatibility.
     */
    static final long serialVersionUID = 1059024301538472131L;
    /**
     * The hydrogen atom.
     */
    public static final Atom H = new Hydrogen();
    /**
     * The nitrogen atom.
     */
    public static final Atom N = new Nitrogen();
    /**
     * The oxygen atom.
     */
    public static final Atom O = new Oxygen();
    /**
     * The carbon atom.
     */
    public static final Atom C = new Carbon();
    /**
     * The sulfur atom.
     */
    public static final Atom S = new Sulfur();
    /**
     * The phosphorus atom.
     */
    public static final Atom P = new Phosphorus();
    /**
     * The helium atom.
     */
    public static final Atom He = new Helium();
    /**
     * The phosphorus atom.
     */
    public static final Atom Li = new Lithium();
    /**
     * The sodium atom.
     */
    public static final Atom Na = new Sodium();
    /**
     * The selenium atom.
     */
    public static final Atom Se = new Selenium();
    /**
     * The Iodine atom.
     */
    public static final Atom I = new Iodine();
    
    /**
     * Returns an array of implemented atoms indicated by their short name.
     * 
     * @param includeSelect if true, the first item is set to '- Select -' 
     * @return an array of implemented atoms
     */
    public static String[] getImplementedAtoms(boolean includeSelect) {
        if (includeSelect) {
            return new String[] {"- Select -", "C", "H", "I", "N", "O", "S", "P", "He", "Li", "Na", "Se"};
        } else {
            return new String[] {"C", "H", "I", "N", "O", "S", "P", "He", "Li", "Na", "Se"};
        }
    }

    /**
     * Returns the atom corresponding to the given short name.
     *
     * @param shortName the short name of the atom
     *
     * @return the atom corresponding to the given short name
     */
    public static Atom getAtom(String shortName) {
        if (shortName.equals("H")) {
            return H;
        } else if (shortName.equals("I")) {
            return I;
        } else if (shortName.equals("N")) {
            return N;
        } else if (shortName.equals("O")) {
            return O;
        } else if (shortName.equals("C")) {
            return C;
        } else if (shortName.equals("S")) {
            return S;
        } else if (shortName.equals("P")) {
            return P;
        } else if (shortName.equals("He")) {
            return He;
        } else if (shortName.equals("Li")) {
            return Li;
        } else if (shortName.equals("Na")) {
            return Na;
        } else if (shortName.equals("Se")) {
            return Se;
        }
        throw new UnsupportedOperationException("Atom " + shortName + " not implemented.");
    }

    /**
     * The monoisotopic mass. Access is faster then querying the isotope map.
     */
    protected Double monoisotopicMass;
    /**
     * Map of the isotope masses relative to the monoisotopic peak (+1 for
     * carbon 13).
     */
    protected HashMap<Integer, Double> isotopeMap;
    /**
     * Map of the isotope representative composition of the stable isotopes.
     */
    protected HashMap<Integer, Double> representativeComposition;
    /**
     * The name of the atom.
     */
    protected String name;
    /**
     * The single letter code of the atom.
     */
    protected String letter;

    /**
     * Returns the monoisotopic mass.
     *
     * @return the monoisotopic mass in Da
     */
    public Double getMonoisotopicMass() {
        return monoisotopicMass;
    }

    /**
     * Returns the name of the atom.
     *
     * @return the name of the atom
     */
    public String getName() {
        return name;
    }

    /**
     * Returns the single letter code of the atom.
     *
     * @return the single letter code of the atom
     */
    public String getLetter() {
        return letter;
    }

    /**
     * returns an unsorted list of isotopes for which a mass is available relative to the
     * monoisotopic peak (+1 for carbon 13).
     *
     * @return a list of isotopes for which a mass is available
     */
    public ArrayList<Integer> getImplementedIsotopes() {
        if (isotopeMap != null) {
            return new ArrayList<Integer>(isotopeMap.keySet());
        }
        return new ArrayList<Integer>();
    }

    /**
     * Returns the mass corresponding to the given isotope number. Null if not
     * found.
     *
     * @param isotopeNumber the isotope number of interest relative to the
     * monoisotopic peak (+1 for carbon 13).
     *
     * @return the corresponding mass
     */
    public Double getIsotopeMass(int isotopeNumber) {
        if (isotopeMap != null) {
            return isotopeMap.get(isotopeNumber);
        }
        return null;
    }

    /**
     * Returns the mass difference between the given isotope and the
     * monoisotopic mass.
     *
     * @param isotopeNumber the isotope number relative to the monoisotopic peak
     * (+1 for carbon 13)
     *
     * @return the mass difference between the given isotope and the
     * monoisotopic mass
     */
    public Double getDifferenceToMonoisotopic(int isotopeNumber) {
        Double isotopeMass = null;
        if (isotopeMap != null) {
            isotopeMass = isotopeMap.get(isotopeNumber);
        }
        if (isotopeMass == null) {
            throw new IllegalArgumentException("No isotope mass found for isotope " + isotopeNumber + " of atom " + name + ".");
        }
        return isotopeMass - getMonoisotopicMass();
    }
}
