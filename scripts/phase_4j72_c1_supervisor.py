"""Linux x86-64 inherited seccomp/ptrace exec gate for C1 software runs.

Only exec syscalls generate seccomp stops, avoiding all-syscall tracing.
Fork/clone descendants inherit the filter; EXITKILL terminates them if the
supervisor dies. Without a tracer SECCOMP_RET_TRACE fails exec closed (ENOSYS).
"""

import ctypes as C
import errno
import os
import signal
from pathlib import Path

from scripts.phase_4j72_c1_fence import (
    append_event,
    command_reason,
    environment_policy,
    executable_identity,
)

libc = C.CDLL(None, use_errno=True)
libc.ptrace.restype = C.c_long


def ptrace(request, pid, address=0, data=0):
    result = libc.ptrace(
        C.c_uint(request), C.c_uint(pid), C.c_void_p(address), C.c_void_p(data)
    )
    if result == -1 and C.get_errno():
        raise OSError(C.get_errno(), "C1 ptrace")
    return result


class Registers(C.Structure):
    _fields_ = [
        (n, C.c_ulonglong)
        for n in [
            "r15",
            "r14",
            "r13",
            "r12",
            "rbp",
            "rbx",
            "r11",
            "r10",
            "r9",
            "r8",
            "rax",
            "rcx",
            "rdx",
            "rsi",
            "rdi",
            "orig_rax",
            "rip",
            "cs",
            "eflags",
            "rsp",
            "ss",
            "fs_base",
            "gs_base",
            "ds",
            "es",
            "fs",
            "gs",
        ]
    ]


class Filter(C.Structure):
    _fields_ = [
        ("code", C.c_ushort),
        ("jt", C.c_ubyte),
        ("jf", C.c_ubyte),
        ("k", C.c_uint),
    ]


class Program(C.Structure):
    _fields_ = [("len", C.c_ushort), ("filter", C.POINTER(Filter))]


def kernel_filter():
    # seccomp_data: nr at 0, arch at 4. Kill unsupported ABI including x32.
    instructions = [
        (0x20, 0, 0, 4),
        (0x15, 1, 0, 0xC000003E),
        (0x06, 0, 0, 0x80000000),
        (0x20, 0, 0, 0),
        (0x45, 0, 1, 0x40000000),
        (0x06, 0, 0, 0x80000000),
        (0x15, 2, 0, 59),
        (0x15, 1, 0, 322),
        (0x06, 0, 0, 0x7FFF0000),
        (0x06, 0, 0, 0x7FF00000),
    ]
    array = (Filter * len(instructions))(*(Filter(*i) for i in instructions))
    program = Program(len(array), array)
    if libc.prctl(38, 1, 0, 0, 0) or libc.prctl(22, 2, C.byref(program), 0, 0):
        raise OSError(C.get_errno(), "C1 seccomp installation failed")


def read_exec(pid, registers):
    with open(f"/proc/{pid}/mem", "rb", buffering=0) as memory:

        def string(address):
            if not address:
                return ""
            data = bytearray()
            # Read page boundaries safely, preserving argv exactly as decoded
            # with surrogateescape; oversized/unreadable requests fail closed.
            while len(data) < 131072:
                chunk = os.pread(
                    memory.fileno(), min(4096 - address % 4096, 4096), address
                )
                if not chunk:
                    raise RuntimeError("unreadable exec string")
                if b"\0" in chunk:
                    return os.fsdecode(bytes(data) + chunk.split(b"\0", 1)[0])
                data.extend(chunk)
                address += len(chunk)
            raise RuntimeError("oversized exec string")

        def vector(address):
            words = []
            for i in range(32768):
                raw = os.pread(memory.fileno(), 8, address + i * 8)
                pointer = int.from_bytes(raw, "little")
                if len(raw) != 8:
                    raise RuntimeError("unreadable exec vector")
                if not pointer:
                    return words
                words.append(string(pointer))
            raise RuntimeError("oversized exec vector")

        if registers.orig_rax == 59:
            path, argv, env = registers.rdi, registers.rsi, registers.rdx
        else:
            # execveat dirfd/empty-path semantics are rejected, never guessed.
            return None, [], {}, "execveat_not_allowed"
        environment = {}
        for word in vector(env):
            key, sep, value = word.partition("=")
            if sep:
                if key in environment:
                    return string(path), vector(argv), {}, "duplicate_environment_key"
                environment[key] = value
        return string(path), vector(argv), environment, None


def supervise(argv, env, policy):
    if os.uname().machine != "x86_64":
        raise RuntimeError("C1 kernel gate supports Linux x86-64 only")
    pid = os.fork()
    if pid == 0:
        try:
            ptrace(0, 0)  # TRACEME
            os.kill(os.getpid(), signal.SIGSTOP)
            kernel_filter()
            os.execve(argv[0], argv, env)
        except BaseException:  # noqa: BLE001 -- failure must stop before exec
            os._exit(97)
    waited, status = os.waitpid(pid, 0)
    if waited != pid or not os.WIFSTOPPED(status):
        raise RuntimeError("C1 child failed before installation")
    # TRACEFORK/VFORK/CLONE/EXEC/SECCOMP plus EXITKILL.
    ptrace(0x4200, pid, data=0x2 | 0x4 | 0x8 | 0x10 | 0x80 | 0x100000)
    ptrace(7, pid)
    live, main_code = {pid}, None
    try:
        while live:
            child, status = os.waitpid(-1, 0x40000000)  # __WALL, all traced tasks
            if os.WIFEXITED(status) or os.WIFSIGNALED(status):
                live.discard(child)
                if child == pid:
                    main_code = os.waitstatus_to_exitcode(status)
                continue
            event, sig = status >> 16, os.WSTOPSIG(status)
            if event in (1, 2, 3):
                new = C.c_ulong()
                ptrace(0x4201, child, data=C.addressof(new))
                live.add(new.value)
            if event == 7:  # PTRACE_EVENT_SECCOMP, before exec
                command, identity, environment = [], None, {}
                try:
                    registers = Registers()
                    ptrace(12, child, data=C.addressof(registers))
                    executable, command, environment, reason = read_exec(
                        child, registers
                    )
                    if not reason:
                        cwd = str(Path(f"/proc/{child}/cwd").resolve())
                        identity = executable_identity(executable, environment, cwd)
                        # Check the actual exec filename AND argv[0]; neither can
                        # disguise a forbidden program or startup bypass.
                        reason = command_reason(command, environment, policy, identity)
                        if not reason:
                            reason = command_reason(
                                [executable, *command[1:]],
                                environment,
                                policy,
                                identity,
                            )
                        # execvpe tries nonexistent PATH candidates before the
                        # actual executable. Preserve ENOENT without allowing
                        # an existing unapproved executable to start.
                        missing = (
                            reason == "executable_not_in_frozen_allowlist"
                            and identity["resolved"] is None
                            and not (Path(cwd) / executable).exists()
                        )
                        if missing:
                            reason = None
                except Exception as exc:  # noqa: BLE001 -- unverifiable means deny
                    reason = "unverifiable_exec:" + type(exc).__name__
                stat = Path(f"/proc/{child}/status").read_text()
                parent = int(
                    next(
                        line.split()[1]
                        for line in stat.splitlines()
                        if line.startswith("PPid:")
                    )
                )
                append_event(
                    policy,
                    {
                        "action": "rejected"
                        if reason
                        else "exec_missing"
                        if identity and identity["resolved"] is None
                        else "exec_allowed",
                        "layer": "kernel_exec",
                        "operation": "execve",
                        "command": command,
                        "executable_identity": identity,
                        "pid": child,
                        "parent_process": parent,
                        "reason": reason,
                        "environment_policy": environment_policy(environment),
                    },
                )
                if reason:
                    # Do not replace or edit argv or execute a substitute.
                    # Kill the stopped child before the requested exec syscall.
                    os.kill(child, signal.SIGKILL)
                    continue
            try:
                ptrace(7, child, data=0 if event or sig == signal.SIGSTOP else sig)
            except OSError as exc:
                if exc.errno != errno.ESRCH:
                    raise
    finally:
        for child in live:
            try:
                os.kill(child, signal.SIGKILL)
            except ProcessLookupError:
                pass
    return main_code
