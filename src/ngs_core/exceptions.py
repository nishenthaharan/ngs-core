"""Domain exceptions for :mod:`ngs_core`."""


class NGSCoreError(Exception):
    """Base exception for errors that should be shown to CLI users."""


class FastqFormatError(NGSCoreError):
    """Raised when an input is not valid four-line FASTQ."""


class PairingError(NGSCoreError):
    """Raised when paired-end inputs are not synchronized."""
