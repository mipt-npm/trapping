plugins {
    kotlin("jvm") version "1.4.10"
    application
}

group = "ru.inr.mass"
version = "1.1.0"

description = "Numass trapping simulation"
application {
    mainClassName = "ru.inr.mass.trapping.MainKt"
}

repositories {
    jcenter()
    mavenCentral()
}

dependencies {
    implementation("org.apache.commons:commons-math3:3.6.1")
    implementation("org.apache.commons:commons-rng-simple:1.1")
    implementation("kscience.plotlykt:plotlykt-core:0.2.0")
    implementation("de.m3y.kformat:kformat:0.7")
    testImplementation("junit:junit:4.12")
}

tasks.compileKotlin {
    kotlinOptions.jvmTarget = "1.8"
}

