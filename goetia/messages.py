#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) Camille Scott, 2019
# File   : messages.py
# License: MIT
# Author : Camille Scott <camille.scott.w@gmail.com>
# Date   : 12.03.2020

from dataclasses import dataclass, field, is_dataclass
from enum import Enum, Flag, IntEnum, auto
from typing import List

from mashumaro import DataClassYAMLMixin, DataClassMessagePackMixin

from goetia.utils import require_kwargs


class MessageType(Enum):
    Interval = auto()
    SampleStarted = auto()
    SampleFinished = auto()
    SampleSaturated = auto()
    Error = auto()
    DistanceCalc = auto()
    EndStream = auto()


class AllMessages:
    pass


@require_kwargs
@dataclass(frozen=True)
class Interval(AllMessages, DataClassYAMLMixin, DataClassMessagePackMixin):
    t: int
    sequence: int
    sample_name: str
    file_names: List[str]

    seconds_elapsed_total: float = 0
    seconds_elapsed_sample: float = 0
    seconds_elapsed_interval: float = 0
    modulus: int = 0
    msg_type: MessageType = MessageType.Interval


@require_kwargs
@dataclass(frozen=True)
class SampleStarted(AllMessages, DataClassYAMLMixin, DataClassMessagePackMixin):
    t: int
    sequence: int
    sample_name: str
    file_names: List[str]

    seconds_elapsed_total: float = 0
    msg_type: MessageType = MessageType.SampleStarted


@require_kwargs
@dataclass(frozen=True)
class SampleFinished(AllMessages, DataClassYAMLMixin, DataClassMessagePackMixin):
    t: int
    sequence: int
    sample_name: str
    file_names: List[str]

    seconds_elapsed_total: float = 0
    seconds_elapsed_sample: float = 0
    msg_type: MessageType = MessageType.SampleFinished


@require_kwargs
@dataclass(frozen=True)
class SampleSaturated(AllMessages, DataClassYAMLMixin, DataClassMessagePackMixin):
    t: int
    sequence: int
    sample_name: str
    file_names: List[str]
    msg_type: MessageType = MessageType.SampleSaturated


@require_kwargs
@dataclass(frozen=True)
class Error(AllMessages, DataClassYAMLMixin, DataClassMessagePackMixin):
    t: int
    sequence: int
    sample_name: str
    error: str
    file_names: List[str]
    msg_type: MessageType = MessageType.Error


@require_kwargs
@dataclass(frozen=True)
class DistanceCalc(AllMessages, DataClassYAMLMixin, DataClassMessagePackMixin):
    t: int
    sequence: int
    sample_name: str
    stat_type: str
    distance: float
    stat: float
    file_names: List[str]
    msg_type: MessageType = MessageType.DistanceCalc


@require_kwargs
@dataclass(frozen=True)
class EndStream(AllMessages, DataClassYAMLMixin, DataClassMessagePackMixin):
    t: int
    sequence: int

    seconds_elapsed_total: float = 0
    msg_type: MessageType = MessageType.EndStream
