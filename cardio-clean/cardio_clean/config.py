# -*- coding: utf-8 -*-

"""
Объект для загрузки конфига
"""

import logging
import os
import yaml


class YAMLConfigReader(object):
    """Object representing YAML settings as a dict-like object
    with values as fields. Settings files are loaded locally from
    config.yaml
    """

    def __init__(self, custom_config_file=None):
        config_files = []  # потенциально может быть несколько файлов
        project_path = os.path.dirname(__file__)
        if not custom_config_file:
            project_config_file = os.path.join(
                project_path,
                "config.yaml"
            )
            config_files.append(project_config_file)
        else:
            config_files.append(custom_config_file)

        self.config = {}
        for sf in config_files:
            try:
                self.update_from_file(sf)
            except Exception as e:
                logging.error(
                    u"Error while reading config file {0}: {1}".format(
                        sf,
                        e
                    )
                )

    def update(self, dct):
        """Update config from dict

        :param dct: dict
        """
        self.config.update(dct)

    def update_from_file(self, path):
        """Update config from YAML file
        """
        with open(path, "r") as custom_config:
            self.config.update(
                yaml.load(custom_config.read())
            )

    def dump(self):
        """Dump config to YAML string
        """
        return yaml.dump(self.config)

    def __getattr__(self, name):
        return self.config.get(name, None)

    def __repr__(self):
        return "<config object>"


config = YAMLConfigReader()
