## Trapping simulation

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
