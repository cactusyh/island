"""Explicit single-thread, non-reentrant OPLS Context ownership."""

from threading import Lock, current_thread

from island.exceptions import (
    EvaluationError,
    EvaluationInputError,
    EvaluationSessionBusyError,
    EvaluationSessionClosedError,
    EvaluationUnavailableError,
)

from .openmm import OpenMMBoundPotential, _openmm, _OpenMMResources
from .oplsaa import OPLSSinglePointEvaluator, _cleanup


class OPLSEvaluationSession(OpenMMBoundPotential):
    """One exclusive Context; fresh verification remains usable after invalidation."""

    def __init__(self, evaluator):
        if type(evaluator) is not OPLSSinglePointEvaluator:
            raise EvaluationInputError("Require OPLSSinglePointEvaluator")
        self._evaluator = evaluator
        self._owner, self._lock = current_thread(), Lock()
        self._closed, self._entered = False, False
        self._resources = None
        mm, _ = _openmm()
        try:
            self._resources = _OpenMMResources(evaluator, mm)
        except Exception as error:
            self._closed = True
            raise EvaluationUnavailableError(
                f"Cannot open OPLS session: {error}"
            ) from error

    @property
    def closed(self):
        return self._closed

    @property
    def model_fingerprint(self):
        return self._evaluator.model_fingerprint

    @property
    def parameter_fingerprint(self):
        return self._evaluator.parameter_fingerprint

    def __copy__(self):
        raise EvaluationInputError("Live sessions cannot be copied")

    def __deepcopy__(self, memo):
        return self.__copy__()

    def _acquire(self):
        if current_thread() is not self._owner:
            raise EvaluationSessionBusyError("Session belongs to its creating thread")
        if not self._lock.acquire(blocking=False):
            raise EvaluationSessionBusyError("Session is in use")

    def _check_open(self):
        if self._closed:
            raise EvaluationSessionClosedError("OPLS session is closed")

    def _dispose(self, original=None):
        self._closed = True
        resources, self._resources = self._resources, None
        if resources is not None:
            _cleanup(resources, original)

    def close(self):
        self._acquire()
        try:
            self._dispose()
        finally:
            self._lock.release()

    def __enter__(self):
        self._acquire()
        try:
            self._check_open()
            if self._entered:
                raise EvaluationSessionBusyError(
                    "Session context manager is not reentrant"
                )
            self._entered = True
            return self
        finally:
            self._lock.release()

    def __exit__(self, _type, value, _traceback):
        self._acquire()
        try:
            self._dispose(value)
        finally:
            self._lock.release()

    def validate_system(self, system):
        self._acquire()
        try:
            self._check_open()
            self._evaluator.validate_system(system)
        finally:
            self._lock.release()

    def evaluate_fresh(self, coordinates=None, *, coordinate_unit="angstrom"):
        self._acquire()
        try:
            return self._evaluator.evaluate_fresh(
                coordinates, coordinate_unit=coordinate_unit
            )
        finally:
            self._lock.release()

    def evaluate(self, coordinates=None, *, coordinate_unit="angstrom"):
        self._acquire()
        try:
            self._check_open()
            xyz = self._evaluator._coordinates(coordinates, coordinate_unit)
            mm, unit = _openmm()
            try:
                return self._evaluator._evaluate_context(
                    xyz, self._resources.context, mm, unit
                )
            except BaseException as error:
                self._dispose(error)
                if not isinstance(error, Exception) or isinstance(
                    error, EvaluationError
                ):
                    raise
                raise EvaluationError(
                    f"OPLS session evaluation failed: {error}"
                ) from error
        finally:
            self._lock.release()
