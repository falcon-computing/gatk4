package com.falconcomputing.genomics;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.net.URL;
import java.util.HashSet;
import java.util.Set;
/**
 * Loads native libraries from the classpath, usually from a jar file.
 */
public final class NativeLibraryLoader {
    private static final Logger logger = LogManager.getLogger(NativeLibraryLoader.class);
    private static final String USE_LIBRARY_PATH = "USE_LIBRARY_PATH";
    private static final Set<String> loadedLibraries = new HashSet<String>();

    // make sure only check license first time library is loaded
    private static boolean licenseChecked  = false; 

    /**
     * Tries to load the native library from the classpath, usually from a jar file. <p>
     *
     * If the USE_LIBRARY_PATH environment variable is defined, the native library will be loaded from the
     * java.library.path instead of the classpath.
     *
     * @param tempDir  directory where the native library is extracted or null to use the system temp directory
     * @param libraryName  name of the shared library without system dependent modifications
     * @return true if the library was loaded successfully, false otherwise
     */
    public static synchronized boolean load(File tempDir, String libraryName) {
        if (loadedLibraries.contains(libraryName)) {
          logger.info(libraryName+" is already loaded");
          return true;
        }
        else if (licenseChecked) {
          // directly return false if license is already checked
          // this implies that library is already loaded
          logger.info(libraryName+" is loaded but unlicensed");
          return false;
        }

        final String systemLibraryName = System.mapLibraryName(libraryName);

        // load from the java.library.path
        // disabled for now
        if (false && System.getenv(USE_LIBRARY_PATH) != null) {
            final String javaLibraryPath = System.getProperty("java.library.path");
            try {
                logger.warn(String.format("OVERRIDE DEFAULT: Loading %s from %s", systemLibraryName, javaLibraryPath));
                logger.warn(String.format("LD_LIBRARY_PATH = %s", System.getenv("LD_LIBRARY_PATH")));
                System.loadLibrary(libraryName);
                return true;
            } catch (Exception|Error e) {
                logger.warn(String.format("OVERRIDE DEFAULT: Unable to load %s from %s", systemLibraryName, javaLibraryPath));
                return false;
            }
        }

        // load from the java classpath
        final String resourcePath = "native/" +  systemLibraryName;
        final URL inputUrl = NativeLibraryLoader.class.getResource(resourcePath);
        if (inputUrl == null) {
            logger.warn("Unable to find native library: " + resourcePath);
            return false;
        }
        logger.info(String.format("Loading %s", systemLibraryName));

        try {
            final File temp = File.createTempFile(FilenameUtils.getBaseName(resourcePath),
                    "." + FilenameUtils.getExtension(resourcePath), tempDir);
            FileUtils.copyURLToFile(inputUrl, temp);
            temp.deleteOnExit();
            logger.debug(String.format("Extracting %s to %s", systemLibraryName, temp.getAbsolutePath()));
            System.load(temp.getAbsolutePath());
        } catch (Exception|Error e) {
            logger.warn(String.format("Unable to load %s from %s (%s)", systemLibraryName, resourcePath, e.getMessage()));
            return false;
        }

        // check if license is valid
        boolean verified = false;
        if (!licenseChecked) {
          if (license_verify()) {
            loadedLibraries.add(libraryName);
            verified = true;
          }
          licenseChecked = true;
        }
        return verified;
    }

    private static native boolean license_verify();
}
