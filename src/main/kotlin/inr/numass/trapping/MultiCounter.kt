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
package inr.numass.trapping

import java.io.PrintStream
import java.lang.Integer.valueOf
import java.util.HashMap

/**
 * TODO есть объект MultiDimensionalCounter, исползовать его?
 *
 * @author Alexander Nozik
 * @version $Id: $Id
 */
class MultiCounter(private var name: String) {

    private var counts = HashMap<String, Int>()

    /**
     *
     * getCount.
     *
     * @param name a [java.lang.String] object.
     * @return a int.
     */
    fun getCount(name: String): Int? {
        return counts[name]
    }

    /**
     *
     * increase.
     *
     * @param name a [java.lang.String] object.
     */
    @Synchronized
    fun increase(name: String) {
        if (counts.containsKey(name)) {
            val count = counts[name]
            counts.remove(name)
            counts[name] = count!! + 1
        } else {
            counts[name] = valueOf(1)
        }

    }

    /**
     *
     * print.
     *
     * @param out a [java.io.PrintWriter] object.
     */
    fun print(out: PrintStream) {
        out.printf("%nValues for counter %s%n%n", this.name)
        for ((keyName, value) in counts) {

            out.printf("%s : %d%n", keyName, value)
        }
    }

    /**
     *
     * reset.
     *
     * @param name a [java.lang.String] object.
     */
    @Synchronized
    fun reset(name: String) {
        if (counts.containsKey(name)) {
            counts.remove(name)
        }
    }

    /**
     *
     * resetAll.
     */
    @Synchronized
    fun resetAll() {
        this.counts = HashMap()
    }
}
