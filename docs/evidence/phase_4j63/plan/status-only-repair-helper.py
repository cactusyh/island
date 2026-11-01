def optimization_budget_exception(exc):
    """Require a specific cap signal; unrelated failures at callback40 stay failures."""
    if getattr(exc, "reason", None) == "J63_optimization_callback_cap":
        return True
    if (
        type(exc).__name__ != "OptimizationConvergenceError"
        or type(exc).__module__ != "psi4.driver.p4util.exceptions"
    ):
        return False
    trace = exc.__traceback__
    while trace is not None:
        frame = trace.tb_frame
        if (
            frame.f_globals.get("__name__") == "psi4.driver.driver"
            and frame.f_code.co_name == "optimize"
        ):
            context = frame.f_locals
            data = context.get("opt_data")
            error = data.get("error") if isinstance(data, dict) else None
            message = error.get("error_message") if isinstance(error, dict) else None
            return message == "Maximum number of steps exceeded: 40." or (
                context.get("n", 0) > 40 and context.get("opt_status") != "FAILED"
            )
        trace = trace.tb_next
    return False
