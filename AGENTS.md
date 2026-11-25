# AGENTS.md

How to work on this repository.

## Ground rules

1. One script at a time.
   Each iteration produces exactly one script (or one focused piece of work).
   Stop there and wait for the owner to review and commit it before starting the
   next piece. Never build ahead of a review.

2. Keep everything in its folder.
   Scripts, data, and results each have a home (see Directory map). Nothing
   lives outside its designated folder.

3. The owner maintains README.md and HISTORY.md. Do not edit them.

4. The agent does not commit. The owner reviews and commits each piece of work;
   the agent leaves its changes in the working tree for that review.

5. Scripts run, and data is fetched and held, on akkad. SSH in with
   `ssh tunnelAkkad`. This repo is the source of truth for code; akkad executes
   it. (Conda env on akkad: `bio`.)

## Directory map

    src/<workstream>/                     analysis scripts, one subfolder per
                                          workstream; each workstream's config
                                          lives next to its script
    results/<workstream>/<script-name>/   each script's outputs, in a subfolder
                                          named after the script
    data/                                 downloaded source data, akkad-only
                                          (gitignored)
    figures/                              repository figures (owner-maintained)

## Notes

- Code is authored from scratch in this repo; do not port or rely on code from
  prior experiments.
