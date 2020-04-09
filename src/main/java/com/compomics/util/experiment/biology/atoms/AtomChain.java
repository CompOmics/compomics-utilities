package com.compomics.util.experiment.biology.atoms;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

/**
 * A chain of atoms.
 *
 * @author Marc Vaudel
 */
public class AtomChain {

    /**
     * Serial number for backward compatibility.
     */
    static final long serialVersionUID = 2222259572093523514L;
    /**
     * The chain of atoms.
     */
    private ArrayList<AtomImpl> atomChain;
    /**
     * The mass of the atom chain.
     */
    private double mass = -1.0;
    /**
     * Cache for the string value.
     */
    private String stringValue = null;

    /**
     * Creates an empty atom chain.
     */
    public AtomChain() {
        atomChain = new ArrayList<>(4);
    }

    /**
     * Returns an atom chain from the input as string. Atoms are represented by
     * their canonical short name, e.g. C for Carbon, Na for Sodium. The
     * occurrence of a given atom is to be written in parentheses, e.g. C(3)PO
     * is parsed as three C's, one P and one O. No negative values are allowed.
     * The isotope is to be written prior to the atom, e.g. 13C(2)18OP is parsed
     * as two 13C atoms, one 18O, and one P.
     *
     * @param atomChainAsString the atomic chain as a string
     *
     * @return the atom chain represented in the given string
     */
    public static AtomChain getAtomChain(String atomChainAsString) {

        AtomChain atomChain = new AtomChain();

        char[] atomChainAsStringCharArray = atomChainAsString.toCharArray();
        for (int i = 0; i < atomChainAsStringCharArray.length; i++) {
            char character = atomChainAsStringCharArray[i];
            if (character != ' ') {
                if (character == '-') {
                    throw new IllegalArgumentException("Negative isotope found in " 
                            + atomChainAsString + ". Please use the atom number, e.g. 13 for 13C.");
                }
                Integer isotopeNumber = 0;
                int charAsInt = Character.getNumericValue(character);
                // Parse isotope number
                while (charAsInt >= 0 && charAsInt <= 9) {
                    isotopeNumber = 10 * isotopeNumber + charAsInt;
                    i++;
                    if (i == atomChainAsStringCharArray.length) {
                        throw new IllegalArgumentException(
                                "Reached the end of the atom chain while parsing "
                                        + "isotope number in " + atomChainAsString + ".");
                    }
                    character = atomChainAsStringCharArray[i];
                    charAsInt = Character.getNumericValue(character);
                }
                // Parse atom name
                StringBuilder atomName = new StringBuilder();
                atomName.append(character);
                Integer occurrence = null;
                if (i + 1 < atomChainAsStringCharArray.length) {
                    char nextCharacter = atomChainAsStringCharArray[i + 1];
                    while (Character.isLowerCase(nextCharacter)) {
                        atomName.append(nextCharacter);
                        i++;
                        if (i + 1 < atomChainAsStringCharArray.length) {
                            nextCharacter = atomChainAsStringCharArray[i + 1];
                        } else {
                            break;
                        }
                    }
                    if (nextCharacter == '(') {
                        // Parse occurrence in parentheses
                        i++;
                        i++;
                        if (i == atomChainAsStringCharArray.length) {
                            throw new IllegalArgumentException(
                                    "Reached the end of the atom chain while "
                                            + "parsing occurrence of " + atomName 
                                            + " in " + atomChainAsString + ".");
                        }
                        character = atomChainAsStringCharArray[i];
                        if (character == '-') {
                            throw new IllegalArgumentException(
                                    "Negative occurrence found for " 
                                            + atomName + " in " + atomChainAsString + ".");
                        }
                        while (character != ')') {
                            charAsInt = Character.getNumericValue(character);
                            if (charAsInt < 0 || charAsInt > 9) {
                                throw new IllegalArgumentException(
                                        "Encountered unexpected character " + character 
                                                + " while parsing occurrence of " 
                                                + atomName + " in " + atomChainAsString + ".");
                            }
                            if (occurrence == null) {
                                occurrence = charAsInt;
                            } else {
                                occurrence = 10 * occurrence + charAsInt;
                            }
                            i++;
                            if (i == atomChainAsStringCharArray.length) {
                                throw new IllegalArgumentException(
                                        "Reached the end of the atom chain while "
                                                + "parsing occurrence of " + atomName 
                                                + " in " + atomChainAsString + ".");
                            }
                            character = atomChainAsStringCharArray[i];
                        }
                    }
                }
                if (occurrence == null) {
                    occurrence = 1;
                }
                Atom atom = Atom.getAtom(atomName.toString());
                AtomImpl atomImpl = new AtomImpl(atom, isotopeNumber);
                if (isotopeNumber != 0) {
                    isotopeNumber = atomImpl.getIsotopeNumber(isotopeNumber);
                    if (isotopeNumber == null) {
                        throw new UnsupportedOperationException(
                                "An error occurred while parsing atom chain " 
                                        + atomChainAsString + "Isotope " 
                                        + isotopeNumber + " not supported for atom " + atom + ".");
                    }
                    atomImpl.setIsotope(isotopeNumber);
                }
                atomChain.append(atomImpl, occurrence);
            }
        }
        return atomChain;
    }

    /**
     * Appends an atom to the chain of atoms.
     *
     * @param atom a new atom
     */
    public void append(AtomImpl atom) {
        atomChain.add(atom);
        stringValue = null;
    }

    /**
     * Appends an atom to the chain of atoms.
     *
     * @param atom a new atom
     * @param occurrence the number of times this atom should be added
     */
    public void append(AtomImpl atom, int occurrence) {
        if (occurrence < 0) {
            throw new IllegalArgumentException("Negative occurrence");
        }
        for (int i = 0; i < occurrence; i++) {
            atomChain.add(atom);
        }
        stringValue = null;
    }

    /**
     * Returns the atom chain as a list of AtomImpl.
     *
     * @return the atom chain as a list of AtomImpl
     */
    public ArrayList<AtomImpl> getAtomChain() {
        return atomChain;
    }

    /**
     * Returns the mass of the atomic chain as sum of the individual atoms.
     *
     * @return the mass of the atomic chain as sum of the individual atoms
     */
    public double getMass() {
        if (mass == -1.0) {
            estimateMass();
        }
        return mass;
    }

    /**
     * Returns the number of atoms in this atom chain.
     *
     * @return the number of atoms in this atom chain
     */
    public int size() {
        return atomChain.size();
    }

    /**
     * Estimates the mass of the atom chain.
     */
    private synchronized void estimateMass() {
        if (mass == -1.0) {
            mass = atomChain.stream()
                    .mapToDouble(AtomImpl::getMass)
                    .sum();
        }
    }

    /**
     * Sets the string value from the stringValue attribute, sets it from the
     * composition if not set.
     *
     * @param includeSpaces boolean indicating whether spaces should be included
     * between atoms
     * @param includeBrackets boolean indicating whether brackets are included
     * for the number of occurrences of each atom
     * @param alwaysIncludeNumber boolean indicating whether the number is shown
     * even when there is only one occurrence of an atom
     * @param isotopeCurlyBrackets if true, the isotopes are indicated as curly
     * brackets after the atom, e.g. carbon 13 is written as C{13}, if false,
     * the isotopes are indicated as a number before the atom, e.g. carbon 13 is
     * written as 13C
     * @param alwaysRecreate if true, the string value will be recreated even if
     * it already exists (use in order to not overwrite the standard format)
     *
     * @return the string value
     */
    public synchronized String getStringValue(boolean includeSpaces, boolean includeBrackets,
            boolean alwaysIncludeNumber, boolean isotopeCurlyBrackets, boolean alwaysRecreate) {

        if (stringValue != null && !alwaysRecreate) {
            return stringValue;
        }

        HashMap<String, Integer> composition = new HashMap<>(atomChain.size());
        HashMap<String, HashMap<Integer, String>> isotopeMap = new HashMap<>(atomChain.size());

        for (AtomImpl atomImpl : atomChain) {
            String atomImplName = atomImpl.toString(isotopeCurlyBrackets);
            Integer occurrence = composition.get(atomImplName);
            if (occurrence == null) {
                occurrence = 0;
            }
            composition.put(atomImplName, occurrence + 1);
            String atomSymbol = atomImpl.getAtomSymbol();
            HashMap<Integer, String> atomIsotopes = isotopeMap.get(atomSymbol);
            if (atomIsotopes == null) {
                atomIsotopes = new HashMap<>(1);
                isotopeMap.put(atomSymbol, atomIsotopes);
            }
            Integer isotope = atomImpl.getIsotope();
            atomIsotopes.put(isotope, atomImplName);
        }

        StringBuilder compositionAsString = new StringBuilder(composition.size());

        ArrayList<String> atomNames = new ArrayList<>(isotopeMap.keySet());
        Collections.sort(atomNames);

        for (String atomLetter : atomNames) {
            HashMap<Integer, String> atomIsotopes = isotopeMap.get(atomLetter);
            ArrayList<Integer> isotopes = new ArrayList<>(atomIsotopes.keySet());
            Collections.sort(isotopes);
            for (Integer isotope : isotopes) {
                String atomName = atomIsotopes.get(isotope);
                if (includeSpaces && compositionAsString.length() > 0) {
                    compositionAsString.append(" ");
                }
                compositionAsString.append(atomName);
                Integer occurrence = composition.get(atomName);
                if (occurrence > 1 || alwaysIncludeNumber) {
                    if (includeBrackets) {
                        compositionAsString.append("(").append(occurrence).append(")");
                    } else {
                        compositionAsString.append(occurrence);
                    }
                }
            }
        }

        String tempStringValue = compositionAsString.toString();

        if (!alwaysRecreate) {
            stringValue = tempStringValue;
        }

        return tempStringValue;
    }

    /**
     * Returns the occurrence of a given atom in the chain.
     *
     * @param atom the atom of interest
     * @param isotope the isotope to look for
     *
     * @return the occurrence of the atom in this atom chain
     */
    public int getOccurrence(Atom atom, Integer isotope) {
        int occurrence = 0;
        AtomImpl atom1 = new AtomImpl(atom, isotope);
        for (AtomImpl atom2 : atomChain) {
            if (atom1.isSameAs(atom2)) {
                occurrence++;
            }
        }
        return occurrence;
    }

    /**
     * Removes all the occurrences of the given atom.
     *
     * @param atom the atom
     * @param isotope the isotope
     */
    public void remove(Atom atom, Integer isotope) {
        ArrayList<AtomImpl> newAtomChain = new ArrayList<>(atomChain.size());
        AtomImpl atom1 = new AtomImpl(atom, isotope);
        for (AtomImpl atom2 : atomChain) {
            if (!atom1.isSameAs(atom2)) {
                newAtomChain.add(atom2);
            }
        }
        atomChain = newAtomChain;
        mass = -1.0;
        stringValue = null;
    }

    /**
     * Sets the occurrence of a given atom.
     *
     * @param atom the atom
     * @param isotope the isotope number
     * @param occurrence the occurrence
     */
    public void setOccurrence(Atom atom, Integer isotope, Integer occurrence) {
        remove(atom, isotope);
        append(new AtomImpl(atom, isotope), occurrence);
        mass = -1.0;
        stringValue = null;
    }

    /**
     * Indicates whether two atom chains are of the same composition by
     * comparing their string and type. An empty chain is considered to be the
     * same composition as a null chain.
     *
     * @param anotherChain another atom chain
     *
     * @return a boolean indicating whether two atom chains are of the same
     * composition
     */
    public boolean isSameCompositionAs(AtomChain anotherChain) {
        if (!atomChain.isEmpty() && (anotherChain == null || anotherChain.getAtomChain().isEmpty())) {
            return false;
        }
        return anotherChain.toString().equals(toString());
    }

    @Override
    public String toString() {
        if (stringValue == null) {
            return getStringValue(false, true, false, false, false);
        }
        return stringValue;
    }

    @Override
    public AtomChain clone() {
        AtomChain result = new AtomChain();
        for (AtomImpl atom : atomChain) {
            result.append(new AtomImpl(atom.getAtomSymbol(), atom.getIsotope()));
        }
        return result;
    }
}
