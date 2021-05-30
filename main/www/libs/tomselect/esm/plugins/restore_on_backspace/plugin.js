/**
* Tom Select v1.7.3
* Licensed under the Apache License, Version 2.0 (the "License");
*/

import TomSelect from '../../tom-select.js';

/**
 * Plugin: "restore_on_backspace" (Tom Select)
 * Copyright (c) contributors
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this
 * file except in compliance with the License. You may obtain a copy of the License at:
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under
 * the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF
 * ANY KIND, either express or implied. See the License for the specific language
 * governing permissions and limitations under the License.
 *
 */
TomSelect.define('restore_on_backspace', function (options) {
  var self = this;

  options.text = options.text || function (option) {
    return option[self.settings.labelField];
  };

  self.on('item_remove', function (value) {
    if (self.control_input.value.trim() === '') {
      var option = self.options[value];

      if (option) {
        self.setTextboxValue(options.text.call(self, option));
      }
    }
  });
});
//# sourceMappingURL=plugin.js.map
