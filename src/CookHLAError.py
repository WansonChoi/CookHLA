#-*- coding: utf-8 -*-

class CookHLAError(Exception):
    """
    Base Error for CookHLA
    """
    pass



class CookHLAImputationError(CookHLAError):
    """
    Error for Beagle 4.1 implementation (Bash Execution).
    """

    def __init__(self, _message):
        print(_message)



class CookHLAInputPreparationError(CookHLAError):
    """
    Error for Beagle 4.1 implementation (Bash Execution).
    """

    def __init__(self, _message):
        print(_message)