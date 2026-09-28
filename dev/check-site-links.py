"""Check local page, fragment, and asset links in a rendered pkgdown site."""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
import sys


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__(convert_charrefs=True)
        self.ids = set()
        self.links = []
        self.feed(text)

    def handle_starttag(self, tag, attributes):
        attrs = dict(attributes)
        self.ids.update(filter(None, (attrs.get("id"), attrs.get("name"))))
        for key in ("href", "src"):
            if attrs.get(key):
                self.links.append(attrs[key])


root = Path(sys.argv[1]).resolve()
pages = {p.resolve(): Page(p.read_text()) for p in root.rglob("*.html")}
if not pages:
    raise SystemExit(f"No HTML pages found under {root}")
errors = []
for path, page in pages.items():
    for link in page.links:
        url = urlsplit(link)
        if url.scheme or url.netloc or link.startswith("/"):
            continue
        target = (path.parent / unquote(url.path)).resolve() if url.path else path
        if target.is_dir():
            target /= "index.html"
        # Release/development navigation may point to the other published site.
        if not target.is_relative_to(root):
            continue
        if not target.exists():
            errors.append(f"{path.relative_to(root)}: missing {link}")
        elif url.fragment and target in pages and unquote(url.fragment) not in pages[target].ids:
            errors.append(f"{path.relative_to(root)}: missing fragment {link}")
if errors:
    raise SystemExit("\n".join(sorted(set(errors))))
print(f"Checked {len(pages)} HTML pages: all local links, fragments, and assets resolve.")
