```@raw html
<!-- Include badges from README -->
<div class="raw-html-chunk">
```

```@eval
# Include just the badges section from README
import Markdown
content = read(joinpath(dirname(@__FILE__), "..", "..", "README.md"), String)
# Extract the badge section (first few lines until the first empty line)
Markdown.parse(content)
```
