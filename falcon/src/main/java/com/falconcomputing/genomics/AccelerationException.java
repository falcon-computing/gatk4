package com.falconcomputing.genomics;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import java.lang.Exception;
import java.lang.String;


/**
 * Falcon custom exception for acceleration procedures
 */
@SuppressWarnings("serial")
public final class AccelerationException extends Exception {
  public AccelerationException() { super(); }
  public AccelerationException(String msg) { super(msg); }
}
