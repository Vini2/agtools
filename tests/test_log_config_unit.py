#!/usr/bin/env python3

import sys

from agtools.log_config import configure_logger
import agtools.log_config as log_config_module


def test_configure_logger_adds_stdout_sink_by_default(monkeypatch):
    calls = []

    monkeypatch.setattr(log_config_module.logger, "remove", lambda: calls.append(("remove",)))
    monkeypatch.setattr(
        log_config_module.logger,
        "add",
        lambda sink, level: calls.append(("add", sink, level)),
    )

    configure_logger()

    assert calls == [("remove",), ("add", sys.stdout, "INFO")]


def test_configure_logger_adds_file_sink_when_requested(monkeypatch, tmp_path):
    calls = []
    log_file = tmp_path / "agtools.log"

    monkeypatch.setattr(log_config_module.logger, "remove", lambda: calls.append(("remove",)))
    monkeypatch.setattr(
        log_config_module.logger,
        "add",
        lambda sink, level: calls.append(("add", sink, level)),
    )

    configure_logger(str(log_file))

    assert calls == [
        ("remove",),
        ("add", str(log_file), "DEBUG"),
    ]


def test_configure_logger_can_add_stderr_sink(monkeypatch):
    calls = []

    monkeypatch.setattr(log_config_module.logger, "remove", lambda: calls.append(("remove",)))
    monkeypatch.setattr(
        log_config_module.logger,
        "add",
        lambda sink, level: calls.append(("add", sink, level)),
    )

    configure_logger(use_stderr=True)

    assert calls == [("remove",), ("add", sys.stderr, "INFO")]
