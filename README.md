[![DOI](https://zenodo.org/badge/261755622.svg)](https://zenodo.org/badge/latestdoi/261755622)

## Trapping simulation

The code for trapping simulation in the [Troitsk nu-mass experiment](http://mass.inr.ru/unu/index_eng.html).

Some design details are discussed in the [Youtube video](https://youtu.be/gG45wzL3gug).

## Authors

* Simulation: Alexander Nozik (INR RAS, MIPT)
* [Initial scattering code](src/scatter/c): Ferenc Glueck and Sebastian Voecking

## Structure

* Electron scattering code in [Scatter.kt](src/main/kotlin/ru/inr/mass/trapping/Scatter.kt).
* Simulation code in [Simulator.kt](src/main/kotlin/ru/inr/mass/trapping/Simulator.kt).

## Dependencies

The simulation geometry relies on [Commons math](https://commons.apache.org/proper/commons-math/) library.

Intermediate pictures created with [Plotly.kt](https://zenodo.org/badge/latestdoi/186020000).

## Building executable

1. Create a fat jar distribution:
```
./gradlew shadowJar
```
The output file is located in `build/libs/trapping-1.1.0-all.jar`

2. Run cross-sections computations
```
java -cp trapping-1.1.0-all.jar ru.inr.mass.trapping.CrosssectionsKt
```

3. Run simulation:
```
./gradlew run   
```
