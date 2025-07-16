import logging
import sys
from pythonjsonlogger import jsonlogger

# --- Custom JSON Formatter ---

class CustomJsonFormatter(jsonlogger.JsonFormatter):
    def add_fields(self, log_record, record, message_dict):
        # Call the parent class method to get the default fields
        super(CustomJsonFormatter, self).add_fields(log_record, record, message_dict)
        
        # Ensure timestamp is in a standard format (ISO 8601)
        if not log_record.get('timestamp'):
            log_record['timestamp'] = record.created
        
        # Rename default fields for clarity and consistency
        if 'message' in log_record:
            log_record['message'] = record.getMessage()
        
        # Add trace_id if it's available in the record
        if 'trace_id' not in log_record and hasattr(record, 'trace_id'):
            log_record['trace_id'] = record.trace_id

# --- Logger Configuration Function ---

_loggers = {}

def get_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """
    Configures and returns a logger with a JSON formatter.
    Ensures that a logger is configured only once.
    """
    if name in _loggers:
        return _loggers[name]

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False # Prevent duplicate logs in parent loggers

    # Use a custom formatter
    formatter = CustomJsonFormatter(
        '%(timestamp)s %(level)s %(name)s %(message)s %(trace_id)s'
    )

    # Use a stream handler to output to stdout
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    _loggers[name] = logger
    return logger

# --- Example Usage (for demonstration) ---

if __name__ == "__main__":
    # Get a logger for a specific module
    log = get_logger("example_module")

    # Standard log
    log.info("This is a standard log message.")

    # Log with extra context (the way we'll pass trace_id)
    log.info("Tool call started.", extra={'trace_id': 'req_xyz789', 'tool_name': 'calculate_logp'})
    log.warning("Something might be wrong.", extra={'trace_id': 'req_xyz789'})
