from __future__ import annotations
from abc import ABC, ABCMeta, abstractmethod
from itertools import dropwhile
from typing import (
    Any,
    Generic,
    Iterable,
    Protocol,
    SupportsIndex,
    Type,
    TypeVar,
    cast,
    overload,
)
import re


T_NMI = TypeVar("T_NMI", bound="IdentMergeableItem")


class IdentMergeableItem(Protocol):
    def ident(self) -> str:
        ...

    def merge(self: T_NMI, other: T_NMI) -> T_NMI:
        ...

    def copy(self: T_NMI) -> T_NMI:
        ...

T = TypeVar("T", bound="UpdateListD")


class UpdateListD(Generic[T_NMI]):
    data: dict[str, T_NMI]
    __slots__ = ("data",)

    def __init__(self, initial: Iterable[T_NMI] = tuple()) -> None:
        self.data = {v.ident(): v for v in initial}

    def __setitem__(self, k: str, v: T_NMI):
        self.data[k] = v

    @overload
    def __getitem__(self, k: str | SupportsIndex) -> T_NMI:
        ...

    @overload
    def __getitem__(self: T, k: slice) -> T:
        ...

    def __contains__(self, kv: str) -> bool:
        return (kv in self.data.keys())

    def __getitem__(self: T, k: str | SupportsIndex | slice) -> T_NMI | T:
        # How easy! It's just a string!
        if isinstance(k, str):
            try:
                return self.data[k]
            except KeyError:
                # A name may have changed, so let's refresh and check before
                # throwing an error.
                self.refreshnames()
                try:
                    return self.data[k]
                except KeyError:
                    raise KeyError(k) from None
        
        m = self.data.values()
        if isinstance(k, SupportsIndex):
            m = iter(m)
            for _ in range(0, k):
                next(m)
            return next(m)
        else:
            if k.step and k.step < 0:
                m = reversed(m)
                step = -k.step
            else:
                m = iter(m)
                if k.step:
                    step = k.step
                else:
                    step = 1
            mi = iter(enumerate(m))
            r = []
            if isinstance(k, slice):
                if isinstance(k.start, SupportsIndex):
                    if int(k.start) < 0:
                        k.start = len(self.data) - int(k.start)
                    mi = dropwhile(lambda x: x[0] < cast(slice, k).start, mi)
                elif k.start is not None:
                    mi = dropwhile(lambda x: x[1].ident() != cast(slice, k).start, mi)
                i = 0
                if isinstance(k.stop, SupportsIndex):
                    while (x := next(mi))[0] < int(k.stop):
                        if i % step == 0:
                            r.append(x[1])
                        i += 1
                elif k.stop is not None:
                    try:
                        while (x := next(mi))[1].ident() != k.stop:
                            if i % step == 0:
                                r.append(x[1])
                            i += 1
                    except StopIteration:
                        raise KeyError(k.stop) from None
                    if i % step == 0:
                        r.append(x[1])
                else:
                    r += list(x[1] for x in mi)
            return self.__class__(r)

    def __delitem__(self, k: str):
        self.data.__delitem__(k)

    def __len__(self):
        return len(self.data)

    def refreshnames(self):
        for k, v in list(self.data.items()):
            if v.ident() != k:
                del self.data[k]
                self.data[v.ident()] = v #FIXME

    def add(self, v: T_NMI):
        k = v.ident()
        if k in self.data:
            self.data[k] = self.data[k].merge(v)
        else:
            self.data[k] = v.copy()

    def __iter__(self):
        return self.data.values().__iter__()

    def __str__(self) -> str:
        return str(self.data.values()).__str__()

    def __repr__(self) -> str:
        return self.__class__.__name__ + "(" + list(self.data.values()).__repr__() + ")"

    def update(self: T, d: T):
        for v in d:
            self.add(v)

    def aslist(self) -> list[T_NMI]:
        return list(self.data.values())

    def asdict(self) -> dict[str, T_NMI]:
        return self.data.copy()

    def copy(self: T) -> T:
        return self.__class__(self.data.copy().values())

    def __or__(self: T, other: T) -> T:
        a = self.copy()
        a.update(other)
        return a

    def __ior__(self: T, other: T) -> T:
        self.refreshnames()
        self.update(other)
        return self

    def __sub__(self: T, other: T) -> T:
        self.refreshnames()
        out = self.copy()
        for v in other:
            if v.name in out:
                # out[v.name].merge(v)  # FIXME
                del(out[v.name])
        return out


    def search(self: T, regex: str, match: bool = False) -> T:
        r = re.compile(regex)
        return self.__class__(v for v in self.data.values() if r.search(v.ident()))


TS = TypeVar("TS", bound="Serializable")


class Serializable(ABC):
    @classmethod
    @abstractmethod
    def _deserialize(cls: Type[TS], input: Any) -> TS:
        ...

    @abstractmethod
    def _serialize(self) -> Any:
        ...

    @classmethod
    def from_yaml(cls: Type[TS], *args, **kwargs) -> TS:
        import yaml

        d = yaml.safe_load(*args, **kwargs)
        return cls._deserialize(d)

    @classmethod
    def from_toml(cls: Type[TS], *args, **kwargs) -> TS:
        import toml

        d = toml.load(*args, **kwargs)
        return cls._deserialize(d)

    @classmethod
    def from_json(cls: Type[TS], *args, **kwargs) -> TS:
        import json

        d = json.load(*args, **kwargs)
        return cls._deserialize(d)

    def to_yaml(self, *args, **kwargs):
        import yaml

        return yaml.safe_dump(self._serialize(), *args, **kwargs)

    def to_json(self, *args, **kwargs):
        import json

        return json.dump(self._serialize(), *args, **kwargs)

    def to_toml(self, *args, **kwargs):
        import toml

        return toml.dump(self._serialize(), *args, **kwargs)
