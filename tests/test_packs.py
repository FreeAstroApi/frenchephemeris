"""The pack download helper (no network, no pack needed)."""

import pytest

from frenchephemeris import packs

pytestmark = pytest.mark.no_packs


def test_sets_reference_published_packs_with_checksums():
    for names in packs.SETS.values():
        for name in names:
            size, sha, _ = packs.PACKS[name]
            assert size > 0 and len(sha) == 64


def test_data_dir_honours_the_environment(monkeypatch, tmp_path):
    monkeypatch.setenv("FRENCHEPHEMERIS_DATA", str(tmp_path))
    assert packs.data_dir() == tmp_path


def test_open_default_explains_how_to_get_packs(monkeypatch, tmp_path):
    import frenchephemeris

    monkeypatch.delenv("FRENCHEPHEMERIS_PACK", raising=False)
    monkeypatch.setenv("FRENCHEPHEMERIS_DATA", str(tmp_path))
    monkeypatch.setattr(frenchephemeris, "pack_dirs", lambda: [tmp_path])
    with pytest.raises(FileNotFoundError, match="frenchephemeris.packs download"):
        frenchephemeris.default_pack_path()
