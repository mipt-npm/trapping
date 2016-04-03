/* 
 * Copyright 2015 Alexander Nozik.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package inr.numass.trapping;

import java.io.PrintStream;
import java.io.PrintWriter;
import static java.lang.Integer.valueOf;
import java.util.HashMap;
import java.util.Map;

/**
 * TODO есть объект MultiDimensionalCounter, исползовать его?
 *
 * @author Alexander Nozik
 * @version $Id: $Id
 */
public class MultiCounter {

    private HashMap<String, Integer> counts = new HashMap<>();
    String name;

    /**
     * <p>Constructor for MultiCounter.</p>
     *
     * @param name a {@link java.lang.String} object.
     */
    public MultiCounter(String name) {
        this.name = name;
    }

    /**
     * <p>getCount.</p>
     *
     * @param name a {@link java.lang.String} object.
     * @return a int.
     */
    public int getCount(String name) {
        if (counts.containsKey(name)) {
            return counts.get(name);
        } else {
            return -1;
        }
    }

    /**
     * <p>increase.</p>
     *
     * @param name a {@link java.lang.String} object.
     */
    public synchronized void increase(String name) {
        if (counts.containsKey(name)) {
            Integer count = counts.get(name);
            counts.remove(name);
            counts.put(name, count + 1);
        } else {
            counts.put(name, valueOf(1));
        }
        
    }

    /**
     * <p>print.</p>
     *
     * @param out a {@link java.io.PrintWriter} object.
     */
    public void print(PrintStream out) {
        out.printf("%nValues for counter %s%n%n", this.name);
        for (Map.Entry<String, Integer> entry : counts.entrySet()) {

            String keyName = entry.getKey();
            Integer value = entry.getValue();
            out.printf("%s : %d%n", keyName, value);
        }
    }

    /**
     * <p>reset.</p>
     *
     * @param name a {@link java.lang.String} object.
     */
    public synchronized void reset(String name) {
        if (counts.containsKey(name)) {
            counts.remove(name);
        }
    }
    /**
     * <p>resetAll.</p>
     */
    public synchronized void resetAll() {
        this.counts = new HashMap<>();
    }
}
