# Issue tracker: GitHub

Issues and specs live in mpaltsai/scRNA-Explorer.
For issue and label operations, use gh with --repo mpaltsai/scRNA-Explorer explicitly.
Verify the target with gh repo view mpaltsai/scRNA-Explorer.

Read issue bodies, labels, and comments before acting.
Create approved tickets using --body-file for multiline content.
Publish blockers first and use native --blocked-by relationships.
Use --parent when creating children of an existing source issue.
If native relationships are unavailable, record blocking issue
references in the body.

PRs as a request surface: no.
