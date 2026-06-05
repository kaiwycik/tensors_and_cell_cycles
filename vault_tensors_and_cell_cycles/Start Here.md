# Start Here — the SCTO note-taking workflow

This vault is built around the **SCTO** method for research note-taking, described
in Ilya Shabanov's article [*A Note-Taking System for Success in Academia*](https://ilyashabanov.substack.com/p/note-taking-system-for-success-in). This note is a short tutorial: read it once, then start taking notes.

**SCTO** is four stages that move an idea from *something you read* to *something
you can publish*:

> **S**ources → **C**ompendiums → **T**houghts → **O**ntologies

In this vault those four stages are the four numbered folders:

| SCTO stage      | Folder here          | What goes in it |
| --------------- | -------------------- | --------------- |
| **S**ources     | `1. Source Notes`    | One note per paper/document you read |
| **C**ompendiums | `2. Notions`         | Facts & concepts in your own words, linked back to their source |
| **T**houghts    | `3. Thoughts`        | Questions and "what-ifs" that pop up while reading |
| **O**ntologies  | `4. Outlines`        | Polished syntheses — drafts of future writing |

![SCTO overview](https://substackcdn.com/image/fetch/$s_!s5NN!,f_auto,q_auto:good,fl_progressive:steep/https%3A%2F%2Fbucketeer-e05bbc84-baa3-437e-9518-adb32be77984.s3.amazonaws.com%2Fpublic%2Fimages%2Faee4dbaf-b466-412c-b3ce-75383d56266e_1306x907.jpeg)

## The workflow, in one pass

**1 — Capture a Source (`1. Source Notes`).**
Every paper you read gets its own source note. Pull it in automatically with `Cmd+P → Zotero Integration: Create Source Note`, or create one by hand. A source note carries its metadata (DOI, rating, date), a short summary, quotable passages, and the PDF itself. Naming follows **`@citekey`** / *FirstAuthor Year*. This is your ground truth — everything downstream links back here.

**2 — Distil into Notions (`2. Notions`).**
A notion is a single finding or concept, written **in your own words**, with a link to the source it came from: `[[@citekey]]`. Keep them small. As a topic grows, split a notion into more specific ones. Collecting every source's findings per concept is how you get an overall view of what you actually know about a subject.

**3 — Capture Thoughts (`3. Thoughts`).**
Loose ideas and questions that surface while reading — *"What if I apply the method from paper X to the problem in paper Y?"* These preserve your *thinking*, not a to-do list. Revisit them; when a question gets answered, archive it rather than delete it.

**4 — Synthesise into Outlines (`4. Outlines`).**
An outline aggregates notions (and, through them, their sources) crossed with your thoughts into a polished, publishable piece — effectively a chapter draft. Because every notion already links to its source, the references compile almost directly into BibTeX.

**Always keep a traceable link back to the primary source.** Notions link to their source note; outlines are built from notions. Use Obsidian's bidirectional `[[links]]` liberally — that web of links is what lets you trace any claim in an outline back to the paper it came from. You may use the **Knowledge Graph** bookmark to visualize that web.

## Daily notes, meetings & to-dos

Idea: **Write in the folder. Dashboard note gathers everything in one place.** The folder is where notes live; the matching `.md` is a read-only view onto them — never edit the dashboard by hand, it rebuilds itself.

| You write in…    | Auto-gathered into… |
| ---------------- | ------------------- |
| `Daily Notes/`   | `Dailies.md`        |
| `Meeting Notes/` | `Meetings.md`       |

**Daily notes.** One note per day to capture progress. Open today's with `Cmd+P → Daily notes: Open today's daily note`, or click a date in the **Calendar**. It lands in `Daily Notes/`, named by date, pre-filled from the daily-note template. **`Dailies.md`** then embeds every daily note, newest first, so you can scroll your whole journal in one place.

**Meeting notes.** One note per meeting. Just create a new note inside the
`Meeting Notes/` folder — the meeting-note template fires automatically, asks for a title, and renames the file `YYYY-MM-DD - Title`. **`Meetings.md`** embeds them all, newest first.

**To-dos.** A to-do is a checkbox line — `- [ ] thing to do` — that you can add in *any* note, most naturally in today's daily note. Use `Cmd+P → Tasks: Create or edit task` to add due dates, priorities, or recurrence. **`TODOs.md`** then collects every open task across the vault into one prioritized, date-aware list — your single to-do hub.

## Supporting notes

- **`Reading List.md`** — papers whose source note still has `status: unread`.
- **`Stats.md`** — Dataview tables describing your corpus (sources by year, notions
  by source, etc.).

Like the dashboards above, these two are read-only views — open them to *see* your corpus, don't edit them by hand.