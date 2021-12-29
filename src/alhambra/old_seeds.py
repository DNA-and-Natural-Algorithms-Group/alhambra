import copy

from .tilestructures import check_edotparen_consistency, check_edotparen_sequence


class seed_base:
    def create_adapter_sequence_diagram(self, adapterdef):
        # At the moment, there are only three potential structures, so I'm just
        # going to put them here.
        ttype = "tile_adapter_5up"
        if "extra" in adapterdef.keys():
            ttype += "_" + adapterdef["extra"]

        import os

        import pkg_resources
        from lxml import etree

        base_svg = etree.parse(
            pkg_resources.resource_stream(
                __name__, os.path.join("seqdiagrambases", ttype + ".svg")
            )
        )

        try:
            s = adapterdef["seqs"][1]
            l = adapterdef["seqs"][0]
        except KeyError:
            s = adapterdef["fullseqs"][1]
            l = adapterdef["fullseqs"][0]
        name = adapterdef["name"]
        if "extra" not in adapterdef.keys():
            strings = [
                (s[0:5] + "--" + s[5:13])[::-1],
                s[13:21] + "--" + s[21:26],
                (l[0:8] + "--" + l[8:16])[::-1],
                l[16:24][::-1],
                l[24:32],
                l[32:40] + "--" + l[40:48],
                adapterdef["ends"][0],
                adapterdef["ends"][1],
                name,
            ]
        elif adapterdef["extra"] == "1h":
            hp = s[0:13]
            s = s[13:]
            strings = [
                hp[0:5] + "-" + hp[5:7] + "-" + hp[7:9],
                (hp[9:11] + "-" + hp[11:13] + "-" + s[0:5] + "--" + s[5:13])[::-1],
                s[13:21] + "--" + s[21:26],
                (l[0:8] + "--" + l[8:16])[::-1],
                l[16:24][::-1],
                l[24:32],
                l[32:40] + "--" + l[40:48],
                adapterdef["ends"][1],
                name,
            ]
        elif adapterdef["extra"] == "2h":
            strings = [
                (s[0:5] + "--" + s[5:13])[::-1],
                s[13:21] + "--" + s[21:26] + "-" + s[26:28] + "-" + s[28:30],
                (s[30:32] + "-" + s[32:34] + "-" + s[34:39])[::-1],
                (l[0:8] + "--" + l[8:16])[::-1],
                l[16:24][::-1],
                l[24:32],
                l[32:40] + "--" + l[40:48],
                adapterdef["ends"][0],
                name,
            ]

        texts = base_svg.findall("//{http://www.w3.org/2000/svg}text")

        for text, string in zip(texts, strings):
            text.text = string

        return base_svg.xpath(
            "/svg:svg/svg:g", namespaces={"svg": "http://www.w3.org/2000/svg"}
        )[0]


class seed_tileadapts(seed_base):
    needspepper = False

    _mimicadapt = {
        "tile_daoe_3up": {
            "strand": 3,
            "ends": slice(1, 3),
            "structure": "tile_adapter_3up",
        },
        "tile_daoe_5up": {
            "strand": 3,
            "ends": slice(1, 3),
            "structure": "tile_adapter_5up",
        },
        "tile_daoe_5up_2h": {
            "strand": 3,
            "ends": slice(1, 3),
            "structure": "tile_adapter_5up_1h",
        },
        "tile_daoe_5up_3b": {
            "strand": 3,
            "ends": slice(1, 3),
            "structure": "tile_adapter_5up_2h",
        },
        "tile_daoe_doublehoriz_35up_2h3h": {
            "strand": 5,
            "ends": slice(2, 4),
            "structure": "tile_adapter_5up_1h",
        },
        "tile_daoe_doublehoriz_35up_1h2i": {
            "strand": 5,
            "ends": slice(2, 4),
            "structure": "tile_adapter_5up",
        },
        "tile_daoe_doublehoriz_35up": {
            "strand": 5,
            "ends": slice(2, 4),
            "structure": "tile_adapter_5up",
        },
        "tile_daoe_doublevert_35up_4h5h": {
            "strand": 5,
            "ends": slice(2, 4),
            "structure": "tile_adapter_5up_2h",
        },
    }

    def create_adapter_sequences(self, tileset):
        tset = copy.deepcopy(tileset)

        for adapter in tset["seed"]["adapters"]:
            # It's not a tile, it's a...
            mimic = tileset.tiles[adapter["tilebase"]]

            structinfo = self._mimicadapt[mimic["structure"]]

            if ("ends" in adapter.keys()) and (
                adapter["ends"] != mimic.ends[structinfo["ends"]]
            ):
                raise ValueError(
                    "adapter {} and tile base {} ends don't match: adapter has {}, tile has {}".format(
                        adapter["name"],
                        mimic.name,
                        adapter["ends"],
                        mimic.ends[structinfo["ends"]],
                    )
                )
            adapter_strand_short = mimic.strands[-1]
            tile_long_strand = mimic.strands[-2]

            adapter["seqs"] = [  # FIXME: better to use short strand here
                tile_long_strand[0:8]
                + self.cores[adapter["loc"] - 1]
                + tile_long_strand[40:48],
                adapter_strand_short,
            ]

            adapter["structure"] = structinfo["structure"]
        return tset

    def check_consistent(self, tileset):
        for adapt in tileset["seed"]["adapters"]:
            mimic = tileset.tiles[adapt["tilebase"]]
            structinfo = self._mimicadapt[mimic["structure"]]
            if adapt.get("structure"):
                assert adapt.get("structure") == structinfo["structure"]
            else:
                adapt["structure"] = structinfo[
                    "structure"
                ]  # FIXME: should be in init code
            if "ends" in adapt.keys():
                assert adapt["ends"] == mimic.ends[structinfo["ends"]]

    edotparens = {
        "tile_adapter_5up": "8(32.8(+5.16)5.",
        "tile_adapter_5up_2h": "8(32.8(+5.16)7(4.7)",
        "tile_adapter_5up_1h": "8(32.8(+7(4.7)16)5.",
    }

    def check_sequence(self, tileset):
        for p in self.edotparens.values():
            check_edotparen_consistency(p)
        self.check_consistent(tileset)
        for adapt in tileset["seed"]["adapters"]:
            if "seqs" not in adapt.keys():
                continue
            seqs = "+".join(adapt["seqs"])
            check_edotparen_sequence(self.edotparens[adapt["structure"]], seqs)
            assert adapt["seqs"][0][8 : (8 + 32)] == self.cores[adapt["loc"] - 1]


class tallrect_tileadapts(seed_tileadapts):
    cores = [
        "CAGGAACGGTACGCCATTAAAGGGATTTTAGA",
        "CTACATTTTGACGCTCACGCTCATGGAAATAC",
        "CCAGCAGAAGATAAAAAATACCGAACGAACCA",
        "GCCGTCAATAGATAATCAACTAATAGATTAGA",
        "ACTTCTGAATAATGGATGATTGTTTGGATTAT",
        "GAAGATGATGAAACAAAATTACCTGAGCAAAA",
        "CATAGGTCTGAGAGACGTGAATTTATCAAAAT",
        "GAAAAAGCCTGTTTAGGGAATCATAATTACTA",
        "ACGCGCCTGTTTATCAGTTCAGCTAATGCAGA",
        "GCTTATCCGGTATTCTAAATCAGATATAGAAG",
        "AACGTCAAAAATGAAAAAACGATTTTTTGTTT",
        "GCAGATAGCCGAACAATTTTTAAGAAAAGTAA",
        "AGACAAAAGGGCGACAGGTTTACCAGCGCCAA",
        "GCGTCAGACTGTAGCGATCAAGTTTGCCTTTA",
        "GTCAGACGATTGGCCTCAGGAGGTTGAGGCAG",
        "TGAAAGTATTAAGAGGCTATTATTCTGAAACA",
    ]


class triangle_side2(seed_tileadapts):
    # Note: this list is currently reversed because
    # the cadnano file for the triangle has 5' ends of
    # the last staple above 3' ends, while our adapter
    # scheme has 3' ends above 5' ends.
    cores = list(
        reversed(
            [
                "AGAGAGTACCTTTAATCCAACAGGTCAGGATT",
                "TAAGAGGAAGCCCGAAATTGCATCAAAAAGAT",
                "CCCCCTCAAATGCTTTTAAATATTCATTGAAT",
                "TACCAGACGACGATAATATCATAACCCTCGTT",
                "CGTTAATAAAACGAACTGGGAAGAAAAATCTA",
                "TAACAAAGCTGCTCATATTACCCAAATCAACG",
                "ATTGTGTCGAAATCCGTATCATCGCCTGATAA",
                "CGAGGGTAGCAACGGCAAAGACAGCATCGGAA",
                "TCTCCAAAAAAAAGGCTTTTTCACGTTGAAAA",
                "TTTCGTCACCAGTACAGTACCGTAACACTGAG",
                "TCAGTACCAGGCGGATATTAGCGGGGTTTTGC",
                "GATACAGGAGTGTACTATACATGGCTTTTGAT",
                "ACCAGAGCCGCCGCCACGCCACCAGAACCACC",
                "GCGTTTGCCATCTTTTCATAGCCCCCTTATTA",
                "ATGAAACCATCGATAGGCCGGAAACGTCACCA",
                "ATCACCGTCACCGACTTCATTAAAGGTGAATT",
            ]
        )
    )


class tallrect_endadapts(seed_base):
    needspepper = True
    cores = tallrect_tileadapts.cores

    def _create_pepper_input_files(self, seeddef, importlist, compstring):
        for adapter in seeddef["adapters"]:
            adaptertype = "tile_adapter_5up"
            e = [[], []]
            e[0].append("origamicore_{0}".format(adapter["loc"]))
            for end in adapter["ends"]:
                if end == "hp":
                    continue  # skip hairpins, etc that aren't designed by stickydesign
                e[0].append("e_" + end.replace("/", "*"))
                if end[-1] == "/":
                    a = "c_" + end[:-1] + "*"
                else:
                    a = "a_" + end
                e[1].append(a)
            s1 = " + ".join(e[0])
            s2 = " + ".join(e[1])
            if "extra" in adapter.keys():
                adaptertype += "_" + adapter["extra"]
            compstring += "component {} = {}: {} -> {}\n".format(
                adapter["name"], adaptertype, s1, s2
            )
            importlist.add(adaptertype)
        return (importlist, compstring)


class bigseed(seed_tileadapts):
    cores = ["A"] * 37


seedtypes = {
    "tallrect_tileadapts": tallrect_tileadapts(),
    "tallrect_endadapts": tallrect_endadapts(),
    "triangle_side2": triangle_side2(),
    "bigseed": bigseed(),
}
seedtypes["longrect"] = seedtypes["tallrect_tileadapts"]
seedtypes["longrect_old"] = seedtypes["tallrect_endadapts"]
