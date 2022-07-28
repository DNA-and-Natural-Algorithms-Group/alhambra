from __future__ import annotations
import base64
import re
from typing import Any, Optional, Union
from typing_extensions import TypeAlias
import xml.etree.ElementTree as ET
import attrs
from attrs import Factory
from abc import ABC, abstractmethod

Number: TypeAlias = Union[int, float]


class DrawingElement(ABC):
    @abstractmethod
    def to_xml(self) -> ET.Element:
        ...


class XMLElement(DrawingElement, ET.Element):
    def to_xml(self) -> ET.Element:
        return self


@attrs.define()
class Group(DrawingElement):
    elements: list[DrawingElement] = Factory(list)
    id: Optional[str] = None

    def to_xml(self) -> ET.Element:
        e = ET.Element("g")
        if self.id is not None:
            e.attrib["id"] = self.id
        for elem in self.elements:
            e.append(elem.to_xml())
        return e

    def append(self, v: DrawingElement) -> None:
        self.elements.append(v)


class Rectangle(DrawingElement):
    _el: ET.Element

    def to_xml(self):
        return self._el

    def __init__(
        self,
        x: Number = 0,
        y: Number = 0,
        width: Number = 1,
        height: Number = 1,
        /,
        **kwargs,
    ):
        e = ET.Element("rect")
        e.attrib["x"] = str(x)
        e.attrib["y"] = str(y)
        e.attrib["width"] = str(width)
        e.attrib["height"] = str(height)
        for k, v in kwargs.items():
            e.attrib[k] = str(v)
        self._el = e


class Text(DrawingElement):
    _el: ET.Element

    def to_xml(self):
        return self._el

    def __init__(
        self, text: str, size: Number = 10, x: Number = 0, y: Number = 0, /, **kwargs
    ):
        e = ET.Element("text")
        e.attrib["x"] = str(x)
        e.attrib["y"] = str(y)
        e.attrib["font-size"] = str(size)
        e.text = text
        for k, v in kwargs.items():
            e.attrib[k.replace("_", "-")] = str(v)
        self._el = e


class Use(DrawingElement):
    _el: ET.Element

    def to_xml(self):
        return self._el

    def __init__(
        self, id_or_link: str | Any, x: Number = 0, y: Number = 0, /, **kwargs
    ):
        e = ET.Element("use")
        e.attrib["x"] = str(x)
        e.attrib["y"] = str(y)
        if not isinstance(id_or_link, str):
            assert hasattr(id_or_link, "id")
            id = str(id_or_link.id)
        else:
            id = id_or_link
        e.attrib["xlink:href"] = "#" + id
        for k, v in kwargs.items():
            e.attrib[k] = str(v)
        self._el = e


@attrs.define()
class Drawing(DrawingElement):
    width: int
    height: int
    defs: list[DrawingElement] = Factory(list)
    elements: list[DrawingElement] = Factory(list)
    viewBox: Optional[tuple[Number, Number, Number, Number]] = None

    def to_xml(self) -> ET.Element:
        e = ET.Element("svg")

        e.attrib["xmlns"] = "http://www.w3.org/2000/svg"
        e.attrib["xmlns:xlink"] = "http://www.w3.org/1999/xlink"
        e.attrib["width"] = str(self.width)
        e.attrib["height"] = str(self.height)

        if self.viewBox is not None:
            e.attrib["viewBox"] = " ".join(str(i) for i in self.viewBox)

        if self.defs:
            s = ET.SubElement(e, "defs")
            for se in self.defs:
                s.append(se.to_xml())
        for se in self.elements:
            e.append(se.to_xml())
        return e

    def save_svg(self, filename: str):
        e = self.to_xml()

        d = ET.ElementTree(e)

        d.write(filename)

    saveSvg = save_svg  # for backwards compatibility

    def to_et(self) -> ET.ElementTree:
        return ET.ElementTree(self.to_xml())

    def to_string(self) -> str:
        e = self.to_xml()
        return ET.tostring(e, encoding="unicode", xml_declaration=True)

    def to_bytes(self) -> bytes:
        e = self.to_xml()
        return ET.tostring(e, encoding="utf-8", xml_declaration=True)

    def _repr_svg_(self):
        return self.to_string()

    # From https://github.com/cduck/drawSvg
    def _repr_html_(self):
        prefix = b"data:image/svg+xml;base64,"
        data = base64.b64encode(self.to_bytes())
        src = (prefix + data).decode(encoding="ascii")
        return '<img src="{}">'.format(src)
