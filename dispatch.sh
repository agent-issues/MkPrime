#!/usr/bin/env bash
# Thin wrapper — delegates to the global dispatch skill.
exec bash "$HOME/.claude/skills/dispatch/dispatch.sh" "$@"
