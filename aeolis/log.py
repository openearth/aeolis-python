import logging

# package modules
from utils import *


logging.basicConfig(filename='aeolis.log',
                    format='%(asctime)-15s %(name)-8s %(levelname)-8s %(message)s')


class Logger():
    '''Logger class

    The logger class support delayed logging of error, warning, info
    and debug messages. Messages logged through this class are
    displayed only upon flush. A user-defined exception can be raised
    if an error is among the flushed messages.

    Examples
    --------
    >>> logger = log.Logger()
    >>> # This error is thrown immediately
    ... log.error('Missing required parameters [%s]' % ', '.join(missing)) 
    >>> # These errors are thrown upon flush
    ... logger.error('Incomplete bathymerty definition')
    ... logger.error('Invalid wind definition file')
    ... logger.error('Wind definition file too short')
    >>> logger.flush(exception=ValueError('Invalid model configuration'))

    '''

    messages = []

    
    def __init__(self):
        pass


    def __getitem__(self, s):
        return self.message[s]
    

    def error(self, msg):
        logging.error(msg)
        return self.add(msg, 'error')


    def warn(self, msg):
        logging.warn(msg)
        return self.add(msg, 'warn')


    def info(self, msg):
        logging.info(msg)
        return self.add(msg, 'info')


    def debug(self, msg):
        logging.debug(msg)
        return self.add(msg, 'debug')


    def add(self, msg, type):
        self.messages.append((msg, type))
        return self


    def flush(self, types=None, exception=None):
        if types is not None:
            types = makeiterable(types)
        error = False
        i = 0
        while i < len(self.messages):
            msg, type = self.messages[i]
            if types is None or type in types:
                print '[%s] %s' % (type.upper(), msg)
                self.messages.pop(i)
                if type == 'error':
                    error = True
            else:
                i += 1
        if error and exception is not None:
            raise exception


def error(msg):
    '''Log error and flush immediately'''
    Logger().error(msg).flush(exception=ValueError(msg))


def warn(msg):
    '''Log warning and flush immediately'''
    Logger().warn(msg).flush()


def info(msg):
    '''Log info message and flush immediately'''
    Logger().info(msg).flush()


def debug(msg):
    '''Log debug message and flush immediately'''
    Logger().debug(msg).flush()
