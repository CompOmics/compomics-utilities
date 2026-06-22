package com.compomics.util.experiment.identification;

import com.compomics.util.db.object.ObjectsDB;
import com.compomics.util.db.object.ObjectsCache;
import com.compomics.util.experiment.identification.matches.PeptideMatch;
import com.compomics.util.experiment.identification.matches.ProteinMatch;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.identification.matches_iterators.PeptideMatchesIterator;
import com.compomics.util.experiment.identification.matches_iterators.ProteinMatchesIterator;
import com.compomics.util.experiment.identification.matches_iterators.PsmIterator;
import com.compomics.util.experiment.identification.matches_iterators.SpectrumMatchesIterator;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.experiment.personalization.ExperimentObject;
import com.compomics.util.waiting.WaitingHandler;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeSet;
import java.util.stream.Collectors;

/**
 * This class interacts with the back-end database to manage identification
 * objects.
 *
 * Interacting with the back-end database might cause InterruptedException.
 * These exceptions are passed as runtime exceptions for methods returning
 * identification objects.
 *
 * @author Marc Vaudel
 * @author Dominik Kopczynski
 */
public class Identification extends ExperimentObject {

    /**
     * Empty default constructor
     */
    public Identification() {
        objectsDB = null;
    }
    /**
     * The keys of the objects in the identification.
     */
    private IdentificationKeys identificationKeys = new IdentificationKeys();
    /**
     * The directory where the database stored.
     */
    private String dbDirectory;
    /**
     * The database which will contain the objects.
     */
    private final ObjectsDB objectsDB;

    /**
     * Constructor.
     *
     * @param objectsDB The object database used to store the objects.
     */
    public Identification(
            ObjectsDB objectsDB
    ) {
        this.objectsDB = objectsDB;
    }

    /**
     * Returns the objects database used in this class.
     *
     * @return the objects database used in this class
     */
    public ObjectsDB getObjectsDB() {
        return objectsDB;
    }

    /**
     * Returns the identification keys.
     *
     * @return the identification keys
     */
    public IdentificationKeys getIdentificationKeys() {
        return identificationKeys;
    }

    /**
     * Sets the identification keys.
     *
     * @param identificationKeys the identification keys
     */
    public void setIdentificationKeys(
            IdentificationKeys identificationKeys
    ) {
        this.identificationKeys = identificationKeys;
    }

    /**
     * Returns a map of the spectrum matches keys indexed by spectrum file name.
     *
     * @return a map of the spectrum matches keys indexed by spectrum file name
     */
    public HashMap<String, HashSet<Long>> getSpectrumIdentification() {
        return identificationKeys.spectrumIdentification;
    }

    /**
     * Returns the spectrum files containing spectrum identifications.
     *
     * @return the spectrum files containing spectrum identifications
     */
    public ArrayList<String> getSpectrumFiles() {
        return new ArrayList<>(identificationKeys.spectrumIdentification.keySet());
    }

    /**
     * Returns the spectrum files containing spectrum identifications.
     *
     * @return the spectrum files containing spectrum identifications
     */
    public ArrayList<String> getOrderedSpectrumFileNames() {
        return getSpectrumFiles();
    }

    /**
     * Returns the spectrum identification keys for the given spectrum file as
     * legacy string spectrum keys.
     *
     * @param spectrumFile the spectrum file
     *
     * @return the spectrum identification keys for the given spectrum file
     */
    public HashSet<String> getSpectrumIdentification(String spectrumFile) {

        HashSet<Long> longKeys = identificationKeys.spectrumIdentification.get(spectrumFile);
        HashSet<String> result = new HashSet<>();

        if (longKeys != null) {

            for (Long key : longKeys) {

                SpectrumMatch spectrumMatch = getSpectrumMatch(key);

                if (spectrumMatch != null) {
                    result.add(Spectrum.getSpectrumKey(spectrumMatch.getSpectrumFile(), spectrumMatch.getSpectrumTitle()));
                }
            }
        }

        return result;
    }

    /**
     * Returns whether the given legacy spectrum key exists.
     *
     * @param spectrumKey the spectrum key
     *
     * @return true if a match exists
     */
    public boolean matchExists(String spectrumKey) {
        return contains(getLongSpectrumKey(spectrumKey));
    }

    /**
     * Returns the assumptions for the given legacy spectrum key.
     *
     * @param spectrumKey the spectrum key
     *
     * @return the assumptions map
     */
    public HashMap<Integer, HashMap<Double, ArrayList<SpectrumIdentificationAssumption>>> getAssumptions(String spectrumKey) {

        SpectrumMatch spectrumMatch = getSpectrumMatch(getLongSpectrumKey(spectrumKey));

        return spectrumMatch == null ? null : spectrumMatch.getAssumptionsMap();
    }

    /**
     * Returns the spectrum match for the given legacy spectrum key.
     *
     * @param spectrumKey the spectrum key
     *
     * @return the spectrum match
     */
    public SpectrumMatch getSpectrumMatch(String spectrumKey) {
        return getSpectrumMatch(getLongSpectrumKey(spectrumKey));
    }

    /**
     * Adds spectrum matches using their own keys.
     *
     * @param spectrumMatches the spectrum matches
     */
    public void addSpectrumMatches(Iterable<SpectrumMatch> spectrumMatches) {

        HashMap<Long, Object> matches = new HashMap<>();

        for (SpectrumMatch spectrumMatch : spectrumMatches) {
            matches.put(spectrumMatch.getKey(), spectrumMatch);
        }

        addSpectrumMatches(matches, null, false);
    }

    /**
     * Returns a PSM iterator for the given spectrum file.
     *
     * @param spectrumFile the spectrum file
     * @param displayProgress display progress
     * @param waitingHandler the waiting handler
     *
     * @return a PSM iterator
     */
    public PsmIterator getPsmIterator(String spectrumFile, boolean displayProgress, WaitingHandler waitingHandler) {

        HashSet<Long> keys = identificationKeys.spectrumIdentification.get(spectrumFile);
        long[] keyArray = keys == null ? new long[0] : keys.stream().mapToLong(Long::longValue).toArray();

        return new PsmIterator(new SpectrumMatchesIterator(keyArray, this, waitingHandler, displayProgress));
    }

    /**
     * Returns a default identification reference.
     *
     * @param projectReference the project reference
     * @param sampleReference the sample reference
     * @param replicateNumber the replicate number
     *
     * @return the default identification reference
     */
    public static String getDefaultReference(String projectReference, String sampleReference, int replicateNumber) {
        return projectReference + "_" + sampleReference + "_" + replicateNumber;
    }

    /**
     * Legacy compatibility method. The current database is provided at
     * construction time.
     *
     * @param dbFolder the database folder
     * @param overwrite overwrite existing database
     * @param objectsCache the objects cache
     */
    public void establishConnection(String dbFolder, boolean overwrite, ObjectsCache objectsCache) {
        if (objectsDB != null && objectsCache != null) {
            objectsDB.setObjectCache(objectsCache);
        }
    }

    /**
     * Resolves a legacy string spectrum key to the current long key.
     *
     * @param spectrumKey the legacy spectrum key
     *
     * @return the long key
     */
    private long getLongSpectrumKey(String spectrumKey) {

        String spectrumFile = Spectrum.getSpectrumFile(spectrumKey);
        String spectrumTitle = Spectrum.getSpectrumTitle(spectrumKey);

        return SpectrumMatch.getKey(spectrumFile, spectrumTitle);
    }

    /**
     * Returns the keys of all the spectrum matches in the db.
     *
     * @return the keys of all the spectrum matches in the db
     */
    public HashSet<Long> getSpectrumIdentificationKeys() {

        return getClassObjects(SpectrumMatch.class);

    }

    /**
     * Returns the number of spectrum identifications.
     *
     * @return the number of spectrum identifications
     */
    public int getSpectrumIdentificationSize() {
        return objectsDB.getNumber(SpectrumMatch.class);
    }

    /**
     * Returns the number of objects of a given class
     *
     * @param className the class name of a given class
     *
     * @return the number of objects
     */
    public int getNumber(
            Class className
    ) {
        return objectsDB.getNumber(className);
    }

    /**
     * Returns an iterator of all objects of a given class
     *
     * @param className the class name of a given class
     * @param filters filters for the class
     *
     * @return the iterator
     */
    /*
    public Iterator<?> getIterator(
            Class className,
            String filters
    ) {
        return objectsDB.getObjectsIterator(className, filters);
    }
    */

    /**
     * Returns the keys of the objects of the given class,
     *
     * @param className the class
     *
     * @return the keys of the objects
     */
    public HashSet<Long> getClassObjects(
            Class className
    ) {
        return objectsDB.getClassObjectIDs(className);
    }

    /**
     * Returns the keys of the objects of the given class,
     *
     * @param className the class
     * @param filters the filters
     *
     * @return the keys of the objects
     */
    public HashSet<Long> getClassObjects(
            Class className,
            String filters
    ) {
        return objectsDB.getClassObjectIDs(className, filters);
    }

    /**
     * Loads all objects of the class in cache.
     *
     * @param className the class name
     * @param waitingHandler the waiting handler allowing displaying progress
     * and canceling the process
     * @param displayProgress boolean indicating whether the progress of this
     * method should be displayed on the waiting handler
     *
     * @throws InterruptedException exception thrown if a threading error occurs
     * while interacting with the database
     */
    public void loadObjects(
            Class className,
            WaitingHandler waitingHandler,
            boolean displayProgress
    ) throws InterruptedException {

        objectsDB.loadObjects(
                className,
                waitingHandler,
                displayProgress
        );

    }

    /**
     * Loads all objects of given keys in cache.
     *
     * @param keyList the list of keys of given objects
     * @param waitingHandler the waiting handler allowing displaying progress
     * and canceling the process
     * @param displayProgress boolean indicating whether the progress of this
     * method should be displayed on the waiting handler
     *
     * @throws InterruptedException exception thrown if a threading error occurs
     * while interacting with the database
     */
    public void loadObjects(
            ArrayList<Long> keyList,
            WaitingHandler waitingHandler,
            boolean displayProgress
    ) throws InterruptedException {

        objectsDB.loadObjects(
                keyList,
                waitingHandler,
                displayProgress
        );

    }

    /**
     * Returns an array of all objects of a given list of keys
     *
     * @param longKey the hash key
     *
     * @return the objects
     */
    public Object retrieveObject(
            long longKey
    ) {

        return objectsDB.retrieveObject(longKey);

    }

    /**
     * Returns the spectrum match with the given key.
     *
     * @param key the key of the match
     *
     * @return the spectrum match with the given key
     */
    public SpectrumMatch getSpectrumMatch(
            long key
    ) {

        return (SpectrumMatch) retrieveObject(key);

    }

    /**
     * Returns the peptide match with the given key.
     *
     * @param key the key of the match
     *
     * @return the peptide match with the given key
     */
    public PeptideMatch getPeptideMatch(
            long key
    ) {

        return (PeptideMatch) retrieveObject(key);

    }

    /**
     * Returns the protein match with the given key.
     *
     * @param key the key of the match
     *
     * @return the protein match with the given key
     */
    public ProteinMatch getProteinMatch(
            long key
    ) {

        return (ProteinMatch) retrieveObject(key);

    }

    /**
     * Returns an array of all objects of a given list of keys
     *
     * @param keyList the key list
     * @param waitingHandler the waiting handler allowing displaying progress
     * and canceling the process
     * @param displayProgress boolean indicating whether the progress of this
     * method should be displayed on the waiting handler
     *
     * @return list of objects
     */
    public ArrayList<Object> retrieveObjects(
            Collection<Long> keyList,
            WaitingHandler waitingHandler,
            boolean displayProgress
    ) {

        return objectsDB.retrieveObjects(
                keyList,
                waitingHandler,
                displayProgress
        );
    }

    /**
     * Returns an array of all objects of a given class
     *
     * @param className the class name
     * @param waitingHandler the waiting handler allowing displaying progress
     * and canceling the process
     * @param displayProgress boolean indicating whether the progress of this
     * method should be displayed on the waiting handler
     *
     * @return list of objects
     */
    public ArrayList<Object> retrieveObjects(
            Class className,
            WaitingHandler waitingHandler,
            boolean displayProgress
    ) {

        return objectsDB.retrieveObjects(
                className,
                waitingHandler,
                displayProgress
        );

    }

    /**
     * Adds an object into the database.
     *
     * @param key the key of the object
     * @param object the object
     */
    public void addObject(
            long key,
            Object object
    ) {

        objectsDB.insertObject(key, object);

    }

    /**
     * Adds a list of objects into the database.
     *
     * @param objects the object
     * @param waitingHandler the waiting handler allowing displaying progress
     * and canceling the process
     * @param displayProgress boolean indicating whether the progress of this
     * method should be displayed on the waiting handler
     */
    public void addObjects(
            HashMap<Long, Object> objects,
            WaitingHandler waitingHandler,
            boolean displayProgress
    ) {

        objectsDB.insertObjects(
                objects,
                waitingHandler,
                displayProgress
        );

    }

    /**
     * Updates an object in the database.
     *
     * @param key the key of the object
     * @param object the object
     */
    public void updateObject(
            long key,
            Object object
    ) {

        objectsDB.updateObject(key, object);

    }

    /**
     * Removes an object from the database.
     *
     * @param key the key of the object
     */
    public void removeObject(
            long key
    ) {

        Object object = objectsDB.retrieveObject(key);

        if (object instanceof ProteinMatch) {

            ProteinMatch proteinMatch = (ProteinMatch) object;

            for (String accession : proteinMatch.getAccessions()) {

                HashSet<Long> proteinKeys = identificationKeys.proteinMap.get(accession);

                if (proteinKeys != null) {

                    proteinKeys.remove(key);

                    if (proteinKeys.isEmpty()) {

                        identificationKeys.proteinMap.remove(accession);

                    }
                }
            }

            identificationKeys.proteinIdentification.remove(key);

        }

        objectsDB.removeObject(key);
    }

    /**
     * Checks if database contains a certain object.
     *
     * @param key the key of the object
     *
     * @return true if database contains a certain object otherwise false
     */
    public boolean contains(
            long key
    ) {

        return objectsDB.inDB(key);

    }

    /**
     * Remove a list of objects from the database.
     *
     * @param keys the list of object keys
     * @param waitingHandler the waiting handler allowing displaying progress
     * and canceling the process
     * @param displayProgress boolean indicating whether the progress of this
     * method should be displayed on the waiting handler
     */
    public void removeObjects(
            ArrayList<Long> keys,
            WaitingHandler waitingHandler,
            boolean displayProgress
    ) {

        objectsDB.removeObjects(
                keys,
                waitingHandler,
                displayProgress
        );

    }

    /**
     * Returns the database directory.
     *
     * @return the database directory
     */
    public String getDatabaseDirectory() {
        return dbDirectory;
    }

    /**
     * Returns a list of the keys of all encountered proteins.
     *
     * @return the corresponding identification results
     */
    public HashSet<Long> getProteinIdentification() {
        return identificationKeys.proteinIdentification;
    }

    /**
     * Returns a list of the keys of all encountered peptides.
     *
     * @return the corresponding identification results
     */
    public HashSet<Long> getPeptideIdentification() {
        return identificationKeys.peptideIdentification;
    }

    /**
     * Adds a peptide match. If an exception occurs when saving to the db it is
     * thrown as runtime exception.
     *
     * @param key the peptide match key
     * @param peptideMatch the peptide match
     */
    public synchronized void addPeptideMatch(
            long key,
            PeptideMatch peptideMatch
    ) {

        identificationKeys.peptideIdentification.add(key);

        objectsDB.insertObject(key, peptideMatch);

    }

    /**
     * Adds a peptide match. If an exception occurs when saving to the db it is
     * thrown as runtime exception.
     *
     * @param peptideMatches the peptide matches in a map
     */
    public synchronized void addPeptideMatches(
            HashMap<Long, Object> peptideMatches
    ) {

        identificationKeys.peptideIdentification.addAll(peptideMatches.keySet());

        objectsDB.insertObjects(peptideMatches, null, false);

    }

    /**
     * Adds a spectrum matches into the database.
     *
     * @param spectrumMatches the spectrum matches
     * @param waitingHandler the waiting handler allowing displaying progress
     * and canceling the process
     * @param displayProgress boolean indicating whether the progress of this
     * method should be displayed on the waiting handler
     */
    public synchronized void addSpectrumMatches(
            HashMap<Long, Object> spectrumMatches,
            WaitingHandler waitingHandler,
            boolean displayProgress
    ) {

        for (Entry<Long, Object> entry : spectrumMatches.entrySet()) {

            long key = entry.getKey();
            SpectrumMatch spectrumMatch = (SpectrumMatch) entry.getValue();

            String fileName = spectrumMatch.getSpectrumFile();

            HashSet<Long> fileKeys = identificationKeys.spectrumIdentification.get(fileName);

            if (fileKeys == null) {

                fileKeys = new HashSet<>();
                identificationKeys.spectrumIdentification.put(fileName, fileKeys);

            }

            fileKeys.add(key);

        }

        objectsDB.insertObjects(
                spectrumMatches,
                waitingHandler,
                displayProgress
        );

    }

    /**
     * Adds a protein match. If an exception occurs when saving to the db it is
     * thrown as runtime exception.
     *
     * @param key The match key.
     * @param proteinMatch The protein match.
     */
    public synchronized void addProteinMatch(
            long key,
            ProteinMatch proteinMatch
    ) {

        fillProteinMap(key, proteinMatch);

        identificationKeys.proteinIdentification.add(key);

        objectsDB.insertObject(key, proteinMatch);

    }

    /**
     * Adds a protein match. If an exception occurs when saving to the db it is
     * thrown as runtime exception.
     *
     * @param proteinMatches Map of the protein matches indexed by key.
     */
    public synchronized void addProteinMatches(
            HashMap<Long, Object> proteinMatches
    ) {

        proteinMatches.entrySet().forEach(
                entry -> fillProteinMap(entry.getKey(), (ProteinMatch) entry.getValue())
        );

        identificationKeys.proteinIdentification.addAll(proteinMatches.keySet());

        objectsDB.insertObjects(proteinMatches, null, false);

    }

    /**
     * Fills the protein map using the given protein match.
     *
     * @param key The match key.
     * @param proteinMatch The protein match.
     */
    private synchronized void fillProteinMap(
            long key,
            ProteinMatch proteinMatch
    ) {

        for (String proteinAccession : proteinMatch.getAccessions()) {

            HashSet<Long> proteinMatchKeys = identificationKeys.proteinMap.get(proteinAccession);

            if (proteinMatchKeys == null) {

                proteinMatchKeys = new HashSet<>(1);
                identificationKeys.proteinMap.put(proteinAccession, proteinMatchKeys);

            }

            proteinMatchKeys.add(key);

        }

    }

    /**
     * Returns a map of all the protein matches which can be ascribed to a
     * protein indexed by its accession.
     *
     * @return a map of all the protein matches which can be ascribed to a
     * protein indexed by its accession.
     */
    public HashMap<String, HashSet<Long>> getProteinMap() {
        return identificationKeys.proteinMap;
    }

    /**
     * Closes the database connection.
     * 
     * @param saveCache save the cache
     */
    public void close(boolean saveCache) {
        objectsDB.close(saveCache);
    }

    /**
     * Closes the database connection and saves the cache.
     */
    public void close() {
        close(true);
    }

    /**
     * Indicates whether the connection to the DB is active.
     *
     * @return true if the connection to the DB is active
     */
    public boolean isConnectionActive() {
        return objectsDB.isConnectionActive();
    }

    /**
     * Returns the keys of the protein matches where a peptide can be found.
     *
     * @param peptideKey the peptide key
     *
     * @return the keys of the protein matches
     */
    public TreeSet<Long> getProteinMatches(
            long peptideKey
    ) {

        return getPeptideMatch(peptideKey)
                .getPeptide()
                .getProteinMapping()
                .navigableKeySet()
                .stream()
                .filter(
                        accession -> identificationKeys.proteinMap.containsKey(accession)
                )
                .flatMap(
                        accession -> identificationKeys.proteinMap.get(accession).stream()
                )
                .filter(
                        proteinKey -> getProteinMatch(proteinKey)
                                .containsPeptide(peptideKey)
                )
                .collect(
                        Collectors.toCollection(TreeSet::new)
                );
    }

    /**
     * Returns a spectrum matches iterator for a given key list.
     *
     * @param spectrumMatches the keys of the spectra to iterate
     * @param waitingHandler the waiting handler
     *
     * @return a spectrum matches iterator
     */
    public SpectrumMatchesIterator getSpectrumMatchesIterator(
            long[] spectrumMatches,
            WaitingHandler waitingHandler
    ) {

        return new SpectrumMatchesIterator(spectrumMatches, this, waitingHandler, false);

    }

    /**
     * Returns a spectrum matches iterator for all SpectrumMatches.
     *
     * @param waitingHandler the waiting handler
     *
     * @return a spectrum matches iterator
     */
    public SpectrumMatchesIterator getSpectrumMatchesIterator(
            WaitingHandler waitingHandler
    ) {

        return new SpectrumMatchesIterator(this, waitingHandler, false);

    }

    /**
     * Returns a psm iterator for all SpectrumMatches.
     *
     * @param waitingHandler the waiting handler
     * @param filters filters for the class
     *
     * @return a peptide matches iterator
     */
    public SpectrumMatchesIterator getSpectrumMatchesIterator(
            WaitingHandler waitingHandler,
            String filters
    ) {

        return new SpectrumMatchesIterator(null, this, waitingHandler, false, filters);

    }

    /**
     * Returns a peptide matches iterator for a given key list.
     *
     * @param peptideKeys the keys of the peptides to iterate
     * @param waitingHandler the waiting handler
     *
     * @return a peptide matches iterator
     */
    public PeptideMatchesIterator getPeptideMatchesIterator(
            long[] peptideKeys,
            WaitingHandler waitingHandler
    ) {

        return new PeptideMatchesIterator(peptideKeys, this, waitingHandler, false);

    }

    /**
     * Returns a peptide matches iterator for all PeptideMatches.
     *
     * @param waitingHandler the waiting handler
     *
     * @return a peptide matches iterator
     */
    public PeptideMatchesIterator getPeptideMatchesIterator(
            WaitingHandler waitingHandler
    ) {

        return new PeptideMatchesIterator(this, waitingHandler, false);

    }

    /**
     * Returns a protein matches iterator for a given key list.
     *
     * @param proteinKeys the keys of the peptides to iterate
     * @param waitingHandler the waiting handler
     *
     * @return a peptide matches iterator
     */
    public ProteinMatchesIterator getProteinMatchesIterator(
            long[] proteinKeys,
            WaitingHandler waitingHandler
    ) {

        return new ProteinMatchesIterator(proteinKeys, this, waitingHandler, false);

    }

    /**
     * Returns a protein matches iterator for all PeptideMatches.
     *
     * @param waitingHandler the waiting handler
     *
     * @return a peptide matches iterator
     */
    public ProteinMatchesIterator getProteinMatchesIterator(
            WaitingHandler waitingHandler
    ) {
        return new ProteinMatchesIterator(this, waitingHandler, false);

    }

    /**
     * Adds a fraction, fractions correspond to the PSM files names. Fractions
     * are ordered alphabetically upon adding of a new fraction.
     *
     * @param fraction the fraction name
     */
    public synchronized void addFraction(
            String fraction
    ) {

        TreeSet orderedFractions = new TreeSet(identificationKeys.fractions);
        orderedFractions.add(fraction);

        setFractions(new ArrayList<>(orderedFractions));

    }

    /**
     * Returns the fractions.
     *
     * @return the fractions
     */
    public ArrayList<String> getFractions() {
        
        return identificationKeys.fractions;
    
    }

    /**
     * Sets the fractions.
     *
     * @param fractions the fractions
     */
    public synchronized void setFractions(
            ArrayList<String> fractions
    ) {

        identificationKeys.fractions = fractions;

    }

}

