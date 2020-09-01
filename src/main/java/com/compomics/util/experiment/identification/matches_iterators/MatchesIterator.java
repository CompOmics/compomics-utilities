package com.compomics.util.experiment.identification.matches_iterators;

import com.compomics.util.db.object.DbObject;
import com.compomics.util.experiment.identification.Identification;
import com.compomics.util.waiting.WaitingHandler;
import java.util.ArrayList;
import java.util.concurrent.Semaphore;
import java.util.Iterator;

/**
 * An abstract iterator class to iterate identification matches. Note: multiple
 * threads are synchronized using semaphores. Should a thread be interrupted,
 * the exception will be sent as RunTimeException. This is because our tools
 * recover from these exceptions similarly as for other unchecked exceptions.
 * Please contact us if you need another/better exception handling.
 *
 * @author Dominik Kopczynski
 * @author Marc Vaudel
 */
public abstract class MatchesIterator {

    /**
     * Empty default constructor
     */
    public MatchesIterator() {
    }
    /**
     * the identification
     */
    private Identification identification = null;
    /**
     * list of potential keys for the iterator, if empty all instances of the
     * class are being iterated
     */
    private long[] keys = null;
    /**
     * current absolute index
     */
    private int index = 0;
    /**
     * absolute number of the objects in the iterator
     */
    private int num = 0;
    /**
     * waiting handler
     */
    private WaitingHandler waitingHandler;
    /**
     * display progress
     */
    private boolean displayProgress;
    /**
     * list of long keys to iterate
     */
    private ArrayList<Long> longKeys = null;
    /**
     * Mutex for the increase of the index.
     */
    private final Semaphore nextMutex = new Semaphore(1);

    /**
     * Constructor.
     *
     * @param className the class name
     * @param identification the identification where to get the matches from
     * @param waitingHandler the waiting handler allowing displaying progress
     * and canceling the process
     * @param displayProgress boolean indicating whether the progress of this
     * method should be displayed on the waiting handler
     */
    public MatchesIterator(Class className, Identification identification, WaitingHandler waitingHandler, boolean displayProgress) {
        this(null, className, identification, waitingHandler, displayProgress, null);
    }

    /**
     * Constructor.
     *
     * @param keys the keys of the objects
     * @param className the className
     * @param identification the identification where to get the matchesloadPs
     * from
     * @param waitingHandler the waiting handler allowing displaying progress
     * and canceling the process
     * @param displayProgress boolean indicating whether the progress of this
     * method should be displayed on the waiting handler
     * @param filters filters for the class
     */
    public MatchesIterator(long[] keys, Class className, Identification identification, WaitingHandler waitingHandler, boolean displayProgress, String filters) {
        
        if (keys != null) {
            
            num = keys.length;
            this.keys = keys;
            
        } else {
            
            longKeys = new ArrayList<>(identification.getClassObjects(className, filters));
             
            
            num = longKeys.size();
            
        }

        index = 0;
        this.identification = identification;
        this.waitingHandler = waitingHandler;
        this.displayProgress = displayProgress;
    }

    /**
     * Returns the next match and updates the buffer. Null if the iterator is
     * done iterating.
     *
     * @return the next match
     */
    public Object nextObject() {

        Object obj = null;
        int currentIndex = getIndex();
        
        if (currentIndex < num) {
            
            if (keys == null) {
                
                obj = identification.retrieveObject(longKeys.get(currentIndex));
                
            } else {
                
                obj = identification.retrieveObject(keys[currentIndex]);
                
            }
        }
        
        return obj;
    }

    /**
     * Returns the index and increases.
     *
     * @return the index
     */
    private int getIndex() {
        
        try {

            nextMutex.acquire();
            int currentIndex = index;
            index++;
            nextMutex.release();
            return currentIndex;

        } catch (InterruptedException e) {

            throw new RuntimeException(e);

        }
    }
}
