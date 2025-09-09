# -------- Colorful printing utility (drop-in, no hard deps) --------
import os
import sys

# ANSI color codes covering common and bright variants
_ANSI_COLOR_CODES = {
    "black": "30", "red": "31", "green": "32", "yellow": "33",
    "blue": "34", "magenta": "35", "cyan": "36", "white": "37",
    "bright_black": "90", "bright_red": "91", "bright_green": "92",
    "bright_yellow": "93", "bright_blue": "94", "bright_magenta": "95",
    "bright_cyan": "96", "bright_white": "97",
}

def _enable_windows_ansi(stream) -> bool:
    """
    Try to enable ANSI escape sequences on Windows consoles.
    Returns True if ANSI is (now or already) supported; False otherwise.
    Strategy:
      1) Try enabling VT100 processing on the console handle via Win32 API.
      2) Fallback: if 'colorama' is installed, initialize it.
    """
    if os.name != "nt":
        return True

    # (1) Attempt to turn on VT processing (Windows 10+)
    try:
        import ctypes
        kernel32 = ctypes.windll.kernel32  # type: ignore[attr-defined]
        # Prefer the stream's handle if possible; otherwise use STD_OUTPUT_HANDLE
        handle = kernel32.GetStdHandle(-11)  # STD_OUTPUT_HANDLE
        mode = ctypes.c_uint32()
        if kernel32.GetConsoleMode(handle, ctypes.byref(mode)):
            ENABLE_VIRTUAL_TERMINAL_PROCESSING = 0x0004
            kernel32.SetConsoleMode(handle, mode.value | ENABLE_VIRTUAL_TERMINAL_PROCESSING)
            return True
    except Exception:
        pass

    # (2) Optional fallback to colorama if available
    try:
        import colorama  # type: ignore
        colorama.just_fix_windows_console()
        return True
    except Exception:
        return False

def _supports_color(stream) -> bool:
    """
    Return True if the given stream is a TTY that likely supports colors.
    - On non-Windows: any TTY is assumed to support ANSI.
    - On Windows: require successful ANSI enablement (see above).
    - If the stream is not a TTY (e.g., piped to file), return False.
    """
    if not hasattr(stream, "isatty") or not stream.isatty():
        return False
    if os.name == "nt":
        return _enable_windows_ansi(stream)
    return True

def print_colorful_message(text: str,
                           color: str = "cyan",
                           bold: bool = False,
                           stream=None) -> None:
    """
    Print colored text to the terminal, gracefully degrading to plain text
    when color is not supported (non-TTY, old terminals, or redirected output).

    Args:
        text:  The message to print.
        color: Color name (e.g., 'blue', 'cyan', 'red', 'bright_blue', ...).
        bold:  If True, apply bold style.
        stream: Target stream (defaults to sys.stdout).

    Behavior:
        - Uses ANSI escape sequences when supported.
        - Unknown color names fall back to 'cyan'.
        - Automatically appends a newline if 'text' does not end with one.
    """
    stream = stream or sys.stdout
    use_color = _supports_color(stream)
    newline = "" if text.endswith("\n") else "\n"

    if use_color:
        code = _ANSI_COLOR_CODES.get(color.lower(), _ANSI_COLOR_CODES["cyan"])
        prefix = f"\033[{'1;' if bold else ''}{code}m"
        suffix = "\033[0m"
        stream.write(f"{prefix}{text}{suffix}{newline}")
    else:
        stream.write(text + newline)
# -------- End of colorful printing utility --------
