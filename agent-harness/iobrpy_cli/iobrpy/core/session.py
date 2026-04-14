"""
session.py - Session management for IOBRpy CLI harness.

This module provides session tracking for interactive CLI sessions,
including history, context, and state persistence.
"""

import json
import shlex
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, List, Any


@dataclass
class CommandHistory:
    """Represents a command in history."""

    command: str
    timestamp: str
    working_dir: str
    exit_code: Optional[int] = None
    output_file: Optional[str] = None
    duration_ms: Optional[int] = None


@dataclass
class SessionContext:
    """Context information for the current session."""

    active_project: Optional[str] = None
    last_command: Optional[str] = None
    last_output_dir: Optional[str] = None
    environment: Dict[str, str] = field(default_factory=dict)
    variables: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> 'SessionContext':
        """Create from dictionary."""
        d = {k: v for k, v in d.items() if v is not None}
        return cls(**d)


class Session:
    """Represents an interactive CLI session."""

    def __init__(self, session_id: str, session_dir: Path):
        """
        Initialize a session.

        Args:
            session_id: Unique session identifier
            session_dir: Directory for session storage
        """
        self.session_id = session_id
        self.session_dir = Path(session_dir) / session_id
        self.session_dir.mkdir(parents=True, exist_ok=True)

        self.context = SessionContext()
        self.history: List[CommandHistory] = []

        self._load()

    def _load(self) -> None:
        """Load session state from disk."""
        context_file = self.session_dir / 'context.json'
        if context_file.exists():
            with open(context_file, 'r') as f:
                self.context = SessionContext.from_dict(json.load(f))

        history_file = self.session_dir / 'history.json'
        if history_file.exists():
            with open(history_file, 'r') as f:
                history_data = json.load(f)
                self.history = [CommandHistory(**h) for h in history_data]

    def _save(self) -> None:
        """Save session state to disk."""
        context_file = self.session_dir / 'context.json'
        with open(context_file, 'w') as f:
            json.dump(self.context.to_dict(), f, indent=2)

        history_file = self.session_dir / 'history.json'
        with open(history_file, 'w') as f:
            json.dump([asdict(h) for h in self.history], f, indent=2)

    def set_active_project(self, project_name: str) -> None:
        """
        Set the active project.

        Args:
            project_name: Name of the active project
        """
        self.context.active_project = project_name
        self._save()

    def get_active_project(self) -> Optional[str]:
        """Get the active project name."""
        return self.context.active_project

    def set_variable(self, name: str, value: Any) -> None:
        """
        Set a session variable.

        Args:
            name: Variable name
            value: Variable value
        """
        self.context.variables[name] = value
        self._save()

    def get_variable(self, name: str, default: Any = None) -> Any:
        """
        Get a session variable.

        Args:
            name: Variable name
            default: Default value if not found

        Returns:
            Variable value or default
        """
        return self.context.variables.get(name, default)

    def record_command(
        self,
        command: str,
        working_dir: str,
        exit_code: Optional[int] = None,
        output_file: Optional[str] = None,
        duration_ms: Optional[int] = None,
    ) -> None:
        """
        Record a command in history.

        Args:
            command: Command string
            working_dir: Working directory
            exit_code: Exit code
            output_file: Output file path
            duration_ms: Duration in milliseconds
        """
        entry = CommandHistory(
            command=command,
            timestamp=datetime.now().isoformat(),
            working_dir=working_dir,
            exit_code=exit_code,
            output_file=output_file,
            duration_ms=duration_ms,
        )
        self.history.append(entry)
        self.context.last_command = command
        self._save()

    def get_history(self, limit: Optional[int] = None) -> List[CommandHistory]:
        """
        Get command history.

        Args:
            limit: Maximum number of entries to return (None for all)

        Returns:
            List of CommandHistory entries
        """
        if limit is None:
            return self.history
        return self.history[-limit:]

    def search_history(self, pattern: str) -> List[CommandHistory]:
        """
        Search command history.

        Args:
            pattern: Search pattern (substring match)

        Returns:
            List of matching CommandHistory entries
        """
        return [
            h for h in self.history
            if pattern.lower() in h.command.lower()
        ]

    def get_last_output_dir(self) -> Optional[str]:
        """Get the last output directory."""
        return self.context.last_output_dir

    def set_last_output_dir(self, path: str) -> None:
        """Set the last output directory."""
        self.context.last_output_dir = path
        self._save()

    def to_dict(self) -> Dict[str, Any]:
        """Convert session to dictionary."""
        return {
            'session_id': self.session_id,
            'path': str(self.session_dir),
            'created_at': datetime.fromtimestamp(self.session_dir.stat().st_ctime).isoformat()
            if self.session_dir.exists() else None,
            'context': self.context.to_dict(),
            'history_count': len(self.history),
        }


class SessionManager:
    """Manages multiple sessions."""

    def __init__(self, sessions_dir: Path):
        """
        Initialize session manager.

        Args:
            sessions_dir: Directory for session storage
        """
        self.sessions_dir = Path(sessions_dir)
        self.sessions_dir.mkdir(parents=True, exist_ok=True)
        self._current_session: Optional[Session] = None

    def create_session(self, session_id: Optional[str] = None) -> Session:
        """
        Create a new session.

        Args:
            session_id: Optional session ID (auto-generated if None)

        Returns:
            New Session instance
        """
        if session_id is None:
            session_id = datetime.now().strftime('session_%Y%m%d_%H%M%S')

        session = Session(session_id, self.sessions_dir)
        self._current_session = session
        return session

    def get_session(self, session_id: str) -> Optional[Session]:
        """
        Get an existing session.

        Args:
            session_id: Session identifier

        Returns:
            Session instance or None if not found
        """
        session_dir = self.sessions_dir / session_id
        if not session_dir.exists():
            return None
        return Session(session_id, self.sessions_dir)

    def get_current_session(self) -> Optional[Session]:
        """Get the current session."""
        return self._current_session

    def set_current_session(self, session: Session) -> None:
        """Set the current session."""
        self._current_session = session

    def list_sessions(self) -> List[Dict[str, Any]]:
        """
        List all sessions.

        Returns:
            List of session info dictionaries
        """
        sessions = []
        for path in self.sessions_dir.iterdir():
            if path.is_dir():
                try:
                    session = Session(path.name, self.sessions_dir)
                    sessions.append(session.to_dict())
                except Exception:
                    pass
        # Sort by creation time, newest first
        sessions.sort(key=lambda x: x.get('created_at', ''), reverse=True)
        return sessions

    def delete_session(self, session_id: str) -> bool:
        """
        Delete a session.

        Args:
            session_id: Session identifier

        Returns:
            True if deleted, False if not found
        """
        session_dir = self.sessions_dir / session_id
        if not session_dir.exists():
            return False

        # Delete session directory
        for item in session_dir.rglob('*'):
            if item.is_file():
                item.unlink()
        for item in sorted(session_dir.rglob('*'), reverse=True):
            if item.is_dir():
                item.rmdir()
        session_dir.rmdir()

        if self._current_session and self._current_session.session_id == session_id:
            self._current_session = None

        return True


def format_command_history(entry: CommandHistory, index: int) -> str:
    """
    Format a command history entry for display.

    Args:
        entry: CommandHistory entry
        index: Entry index

    Returns:
        Formatted string
    """
    status = "OK" if entry.exit_code == 0 else "FAIL" if entry.exit_code is not None else " "
    output_info = f" -> {Path(entry.output_file).name}" if entry.output_file else ""
    duration = f" ({entry.duration_ms}ms)" if entry.duration_ms else ""

    return (
        f"  [{index}] {entry.command}{output_info}{duration}\n"
        f"       {entry.timestamp} | Status: {status}"
    )
