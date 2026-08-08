"""Make sure a pandoc-built .docx declares a content type for every image it embeds.

A .docx is a zip, and Word will only render a part whose type is declared in `[Content_Types].xml`
— either by a `<Default Extension="png" .../>` for the extension or by an `<Override>` for that one
part. Pandoc writes the media files and, depending on the version and on what the source contained,
sometimes writes no `Default` for `png`. Word then shows every figure as an empty frame with a red
cross, which looks like a broken build and is really one missing line of XML.

This adds the missing `Default` entries in place, and says how many images it found either way, so a
figure that silently failed to embed is visible from the build log rather than from opening the file.

Run:  python fix_docx_images.py build/zombi2-manual.docx
"""

from __future__ import annotations

import pathlib
import re
import shutil
import sys
import tempfile
import zipfile

#: extension -> the MIME type Word expects for it
_TYPES = {
    "png": "image/png",
    "jpeg": "image/jpeg",
    "jpg": "image/jpeg",
    "gif": "image/gif",
    "emf": "image/x-emf",
    "wmf": "image/x-wmf",
}


def fix(path: pathlib.Path) -> int:
    """Add any missing image `Default` to ``path``'s `[Content_Types].xml`; return how many."""
    with zipfile.ZipFile(path) as zf:
        names = zf.namelist()
        content_types = zf.read("[Content_Types].xml").decode("utf-8")
        payload = {n: zf.read(n) for n in names}

    media = [n for n in names if n.startswith("word/media/")]
    have = set(re.findall(r'<Default\s+Extension="([^"]+)"', content_types))
    want = {p.rsplit(".", 1)[-1].lower() for p in media} & set(_TYPES)
    missing = sorted(want - have)

    if missing:
        added = "".join(f'<Default Extension="{e}" ContentType="{_TYPES[e]}"/>' for e in missing)
        content_types = content_types.replace("<Types ", "<Types ", 1)
        content_types = re.sub(r"(<Types[^>]*>)", r"\1" + added, content_types, count=1)
        payload["[Content_Types].xml"] = content_types.encode("utf-8")
        with tempfile.NamedTemporaryFile(delete=False, suffix=".docx") as tmp:
            with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as out:
                for name in names:
                    out.writestr(name, payload[name])
        shutil.move(tmp.name, path)

    print(f"  {path.name}: {len(media)} images, "
          + (f"declared {', '.join(missing)}" if missing else "content types already complete"))
    return len(missing)


if __name__ == "__main__":
    fix(pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "build/zombi2-manual.docx"))
