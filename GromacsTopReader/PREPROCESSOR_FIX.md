# Preprocessor Directive Handling

## Issue

The Python parser (`gromacs_to_psf.py`) was encountering `ValueError` when parsing GROMACS topology files containing preprocessor directives:

```
ValueError: invalid literal for int() with base 10: '#ifdef'
```

This occurred in the GO_8AXD files which contain conditional compilation directives:
- `#ifdef` / `#ifndef`
- `#define`
- `#else`
- `#endif`

## Root Cause

The Python parser was only handling `#include` directives and attempting to parse all other lines (including preprocessor directives) as topology data within sections like `[ bonds ]`, `[ angles ]`, etc.

## Solution

Updated both `parse_itp_file()` and `parse_top_file()` to skip all preprocessor directives except `#include`:

```python
# Skip other preprocessor directives (#ifdef, #ifndef, #else, #endif, #define, etc.)
if original_line.strip().startswith('#'):
    continue
```

**Implementation details:**
- Check is performed on the `original_line` (before comment stripping)
- Occurs after `#include` handling but before section header parsing
- Ensures preprocessor directives don't interfere with topology parsing

## C Implementation Behavior

The C implementation (`grotopplugin.c`) already handled this correctly:

1. Checks for `#include` and processes it
2. Strips comments from the line
3. Checks if line is empty or a section header
4. Silently ignores any other lines (including preprocessor directives)

The preprocessor directives don't match the section header pattern `[...]`, so they're naturally skipped in the main parsing loop.

## Test Results

### GO_8AXD System (minimization/system.top)

**Before fix:** `ValueError: invalid literal for int() with base 10: '#ifdef'`

**After fix:**

| Implementation | Atoms | Bonds | Angles | Dihedrals | Impropers | Status |
|----------------|-------|-------|--------|-----------|-----------|--------|
| C | 35,976 | 13,546 | - | - | - | ✅ |
| Python | 35,976 | 13,546 | 11,661 | 2,764 | 450 | ✅ |

Both implementations now successfully parse the file and produce identical atom and bond counts.

## Preprocessor Directives in GROMACS

GROMACS topology files commonly use preprocessor directives for:

1. **Conditional compilation:**
   ```gromacs
   #ifndef FLEXIBLE
   [ constraints ]
   ...
   #else
   [ bonds ]
   ...
   #endif
   ```

2. **Define constants:**
   ```gromacs
   #ifndef POSRES_FC
   #define POSRES_FC 1000
   #endif
   ```

3. **Include files:**
   ```gromacs
   #include "forcefield.itp"
   #include "molecules.itp"
   ```

## Implementation Status

✅ **C Implementation** - Already handled correctly (silently ignores non-include preprocessor directives)

✅ **Python Implementation** - Now aligned with C behavior (explicitly skips non-include preprocessor directives)

## Files Modified

- `gromacs_to_psf.py` (lines 230-232, 320-322)

## Future Considerations

Currently, both implementations treat preprocessor directives as simple text to skip. A more sophisticated implementation could:

1. Evaluate `#define` statements and substitute values
2. Process `#ifdef`/`#ifndef` conditionals
3. Support conditional topology generation

However, for most use cases, simply skipping these directives works well since:
- Most preprocessing is done by `grompp` before our parsers run
- VMD typically loads pre-processed structures
- The goal is to read existing topology files, not compile them
