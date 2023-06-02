.. highlight:: none

.. _topics-commands:

====================================
Command line tool
====================================

Scrapy is controlled through the ``scrapy`` command-line tool, to be referred
here as the "Scrapy tool" to differentiate it from the sub-commands, which we
just call "commands" or "Scrapy commands".

The Scrapy tool provides several commands, for multiple purposes, and each one
accepts a different set of arguments and options.


Configuration settings
====================================

Scrapy will look for configuration parameters in ini-style ``scrapy.cfg`` files
in standard locations:

1. ``/etc/scrapy.cfg`` or ``c:\scrapy\scrapy.cfg`` (system-wide),
2. ``~/.config/scrapy.cfg`` (``$XDG_CONFIG_HOME``) and ``~/.scrapy.cfg`` (``$HOME``)
   for global (user-wide) settings, and
3. ``scrapy.cfg`` inside a Scrapy project's root (see next section).

Settings from these files are merged in the listed order of preference:
user-defined values have higher priority than system-wide defaults
and project-wide settings will override all others, when defined.


Using the ``bust`` tool
====================================

You can start by running the Scrapy tool with no arguments and it will print
some usage help and the available commands::

    Scrapy X.Y - no active project

    Usage:
      scrapy <command> [options] [args]

    Available commands:
      crawl         Run a spider
      fetch         Fetch a URL using the Scrapy downloader
    [...]

The first line will print the currently active project if you're


Available tool commands
====================================

This section contains a list of the available built-in commands with a
description and some usage examples. Remember, you can always get more info
about each command by running::

    bust <command> -h

And you can see all available commands with::

    bust -h



.. Global commands:

.. * :command:`redock`
.. * :command:`dock`
.. * :command:`mol`
.. * :command:`table`
.. * :command:`version`

.. .. command:: redock

.. redock
.. ------------------------------------



.. .. command:: dock

.. dock
.. ------------------------------------


.. .. command:: mol

.. mol
.. ------------------------------------



.. .. command:: table

.. table
.. ------------------------------------


.. .. command:: version

.. version
.. -------

.. * Syntax: ``scrapy version [-v]``
.. * Requires project: *no*

.. Prints the Scrapy version. If used with ``-v`` it also prints Python, Twisted
.. and Platform info, which is useful for bug reports.
