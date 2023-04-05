"""
This file is part of AeoLiS test suit. It uses the Pytest framework.

To run all tests, use:
    pytest

To run only the test in this file, use:
    pytest aeolis/tests/test_utils.py

To run all test in a group (class), use:
    pytest aeolis/tests/test_utils.py::<TestClass>

To run an specific test, use:
    - If in a group:
        pytest aeolis/tests/test_utils.py::<TestClass>::<test-method>
    
    - If not in any group:
        pytest aeolis/tests/test_utils.py::<test-function>
"""

import logging
from aeolis.model import StreamFormatter

class TestStreamFormatter:
    """Test the stream formatter for console output"""

    def test_stream_formatter_parent(self):
        """Test if the stream formatter inherits from the logging.Formatter class"""
        stream_formatter = StreamFormatter()
        assert isinstance(stream_formatter, logging.Formatter)
    
    def test_stream_formatter_info_level(self):
        """Test if the stream formatter returns log message for the info level"""
        logger = logging.getLogger("Test")
        info_record = logger.makeRecord("Test", logging.INFO, 
                                        "Test", 1, "This is a message for INFO level",
                                        None, None)
        stream_formatter = StreamFormatter()
        info_message = stream_formatter.format(info_record)
        
        warning_record = logger.makeRecord("Test warning", logging.WARNING,
                                        "Test", 2, "This is a message for WARNING level",
                                        None, None)

        warning_message = stream_formatter.format(warning_record)

        assert info_message == "This is a message for INFO level"
        assert warning_message == "WARNING: This is a message for WARNING level"

